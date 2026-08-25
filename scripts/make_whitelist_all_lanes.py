import argparse
import gzip
from collections import Counter
from pathlib import Path


DEFAULT_MATRIX_FASTQ = (
    "/data/5th_data/matrix_121324_p4_2nM/test_S1_R1_001.fastq.gz"
)

DEFAULT_TILE_FILES = {
    "1000_lane1": "/data/1000_lane1_xynum.txt",
    "1000_lane2": "/data/1000_lane2_xynum.txt",
    "2000_lane1": "/data/2000_lane1_xynum.txt",
    "2000_lane2": "/data/2000_lane2_xynum.txt",
}

DEFAULT_OUTPUTS = {
    "1000_lane1": "/data/whitelist_1000_lane1_RC.txt",
    "1000_lane2": "/data/whitelist_1000_lane2_RC.txt",
    "2000_lane1": "/data/whitelist_2000_lane1_RC.txt",
    "2000_lane2": "/data/whitelist_2000_lane2_RC.txt",
}


def reverse_complement(seq: str) -> str:
    complement = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )
    return seq.translate(complement)[::-1]


def load_tiles(tile_file: str) -> set[str]:
    tiles: set[str] = set()

    with open(tile_file, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split()

            if parts and parts[0].isdigit():
                tiles.add(parts[0])

    return tiles


def build_tile_to_group(
    group_to_tilefile: dict[str, str]
) -> tuple[dict[str, str], dict[str, set[str]]]:

    group_tiles: dict[str, set[str]] = {}
    tile_to_group: dict[str, str] = {}

    for group, tile_file in group_to_tilefile.items():
        tiles = load_tiles(tile_file)

        if not tiles:
            raise ValueError(
                f"No tiles found in {tile_file} for group {group}"
            )

        group_tiles[group] = tiles

        for tile in tiles:
            if tile in tile_to_group:
                raise ValueError(
                    f"Tile {tile} appears in both "
                    f"{tile_to_group[tile]} and {group}. "
                    "Tile lists must not overlap."
                )

            tile_to_group[tile] = group

    return tile_to_group, group_tiles


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Create 4 lane-specific whitelist files in one pass "
            "over the matrix FASTQ. "
            "Extract bases [start:end), reverse-complement, "
            "and dispatch reads according to the tile ID."
        )
    )

    p.add_argument(
        "--matrix-fastq",
        default=DEFAULT_MATRIX_FASTQ
    )

    p.add_argument(
        "--tile-1000-lane1",
        default=DEFAULT_TILE_FILES["1000_lane1"]
    )

    p.add_argument(
        "--tile-1000-lane2",
        default=DEFAULT_TILE_FILES["1000_lane2"]
    )

    p.add_argument(
        "--tile-2000-lane1",
        default=DEFAULT_TILE_FILES["2000_lane1"]
    )

    p.add_argument(
        "--tile-2000-lane2",
        default=DEFAULT_TILE_FILES["2000_lane2"]
    )

    p.add_argument(
        "--out-1000-lane1",
        default=DEFAULT_OUTPUTS["1000_lane1"]
    )

    p.add_argument(
        "--out-1000-lane2",
        default=DEFAULT_OUTPUTS["1000_lane2"]
    )

    p.add_argument(
        "--out-2000-lane1",
        default=DEFAULT_OUTPUTS["2000_lane1"]
    )

    p.add_argument(
        "--out-2000-lane2",
        default=DEFAULT_OUTPUTS["2000_lane2"]
    )

    p.add_argument(
        "--start-pos",
        type=int,
        default=4,
        help="0-based inclusive start position"
    )

    p.add_argument(
        "--end-pos",
        type=int,
        default=36,
        help="0-based exclusive end position"
    )

    p.add_argument(
        "--orientation",
        choices=["RC", "FWD"],
        default="RC",
        help="Whether to write reverse-complemented or forward barcodes"
    )

    p.add_argument(
        "--summary-tsv",
        default="/data/whitelist_lane_summary.tsv",
        help="Summary TSV with tile counts and matched read counts"
    )

    p.add_argument(
        "--progress-every",
        type=int,
        default=1_000_000,
        help="Print progress every N reads"
    )

    return p.parse_args()


def main() -> None:
    args = parse_args()

    group_to_tilefile = {
        "1000_lane1": args.tile_1000_lane1,
        "1000_lane2": args.tile_1000_lane2,
        "2000_lane1": args.tile_2000_lane1,
        "2000_lane2": args.tile_2000_lane2,
    }

    group_to_output = {
        "1000_lane1": args.out_1000_lane1,
        "1000_lane2": args.out_1000_lane2,
        "2000_lane1": args.out_2000_lane1,
        "2000_lane2": args.out_2000_lane2,
    }

    tile_to_group, group_tiles = build_tile_to_group(
        group_to_tilefile
    )

    for group, outfile in group_to_output.items():
        Path(outfile).parent.mkdir(
            parents=True,
            exist_ok=True
        )

    print("=== lane-specific whitelist creation start ===")
    print(f"Matrix FASTQ: {args.matrix_fastq}")
    print(
        f"Barcode slice: "
        f"[{args.start_pos}:{args.end_pos}] "
        f"({args.end_pos - args.start_pos} bp)"
    )
    print(f"Orientation: {args.orientation}")

    for group in [
        "1000_lane1",
        "1000_lane2",
        "2000_lane1",
        "2000_lane2",
    ]:
        print(
            f"  {group}: "
            f"{len(group_tiles[group])} tiles "
            f"-> {group_to_output[group]}"
        )

    total_reads = 0
    malformed_headers = 0
    unmatched_tiles = 0
    written_counts = Counter()

    writers = {
        group: open(
            outfile,
            "w",
            encoding="utf-8"
        )
        for group, outfile in group_to_output.items()
    }

    try:
        with gzip.open(
            args.matrix_fastq,
            "rt",
            encoding="utf-8",
            errors="replace"
        ) as f_in:

            while True:
                header = f_in.readline().strip()

                if not header:
                    break

                seq = f_in.readline().strip()
                f_in.readline()
                f_in.readline()

                total_reads += 1

                try:
                    tile_id = header.split(":")[4]
                except IndexError:
                    malformed_headers += 1
                    continue

                group = tile_to_group.get(tile_id)

                if group is None:
                    unmatched_tiles += 1
                    continue

                barcode = seq[
                    args.start_pos:args.end_pos
                ]

                if len(barcode) != (
                    args.end_pos - args.start_pos
                ):
                    unmatched_tiles += 1
                    continue

                if args.orientation == "RC":
                    barcode = reverse_complement(
                        barcode
                    )

                writers[group].write(
                    barcode + "\n"
                )

                written_counts[group] += 1

                if (
                    total_reads
                    % args.progress_every
                    == 0
                ):
                    print(
                        f"Processed "
                        f"{total_reads:,} reads | "
                        + ", ".join(
                            f"{g}="
                            f"{written_counts[g]:,}"
                            for g in group_to_output
                        )
                    )

    finally:
        for fh in writers.values():
            fh.close()

    with open(
        args.summary_tsv,
        "w",
        encoding="utf-8"
    ) as out:

        out.write(
            "group\t"
            "tile_file\t"
            "n_tiles\t"
            "output\t"
            "matched_reads\n"
        )

        for group in [
            "1000_lane1",
            "1000_lane2",
            "2000_lane1",
            "2000_lane2",
        ]:
            out.write(
                f"{group}\t"
                f"{group_to_tilefile[group]}\t"
                f"{len(group_tiles[group])}\t"
                f"{group_to_output[group]}\t"
                f"{written_counts[group]}\n"
            )

        out.write(
            f"TOTAL\t-\t-\t-\t"
            f"{sum(written_counts.values())}\n"
        )

        out.write(
            f"UNMATCHED_TILES\t-\t-\t-\t"
            f"{unmatched_tiles}\n"
        )

        out.write(
            f"MALFORMED_HEADERS\t-\t-\t-\t"
            f"{malformed_headers}\n"
        )

        out.write(
            f"TOTAL_READS\t-\t-\t-\t"
            f"{total_reads}\n"
        )

    print("\n=== done ===")

    print(
        f"Total reads:       "
        f" {total_reads:,}"
    )

    print(
        f"Malformed headers: "
        f" {malformed_headers:,}"
    )

    print(
        f"Unmatched tiles:   "
        f" {unmatched_tiles:,}"
    )

    for group in [
        "1000_lane1",
        "1000_lane2",
        "2000_lane1",
        "2000_lane2",
    ]:
        print(
            f"{group:>10}: "
            f"{written_counts[group]:,} "
            f"-> {group_to_output[group]}"
        )

    print(
        f"Summary: "
        f"{args.summary_tsv}"
    )


if __name__ == "__main__":
    main()

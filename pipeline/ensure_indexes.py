"""Add the scalar indexes the paper search uses to filter before ranking.

Does not build a vector index. Search stays an exact comparison of the rows
that pass the filter.
"""

from __future__ import annotations

from pipeline.retrieve import DB_DIR


def main() -> None:
    import lancedb

    database = lancedb.connect(str(DB_DIR))
    table = database.open_table("papers")
    wanted = [
        ("gene_symbols", "LABEL_LIST"),
        ("year", "BTREE"),
        ("month", "BTREE"),
        ("day", "BTREE"),
        ("before_cutoff", "BITMAP"),
    ]
    have = set()
    for index in table.list_indices():
        columns = list(getattr(index, "columns", []) or [])
        kind = str(getattr(index, "index_type", "") or "")
        for column in columns:
            have.add((column, kind.upper()))
        print(f"HAVE {kind} {columns}", flush=True)
    for column, kind in wanted:
        if any(name == column and kind in stored for name, stored in have):
            print(f"SKIP {kind} {column}", flush=True)
            continue
        print(f"BUILD {kind} {column}", flush=True)
        table.create_scalar_index(column, index_type=kind, replace=True)
        print(f"BUILT {kind} {column}", flush=True)
    print("INDEXES_DONE", flush=True)


if __name__ == "__main__":
    main()

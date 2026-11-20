"""Energy import and deterministic conformer selection."""
from .io import finite, read_csv, write_csv


def load_energies(path, catalyst, facet="", aliases=None):
    rows = read_csv(path, ("catalyst", "species", "free_energy"))
    aliases = aliases or {}
    selected, skipped = {}, 0
    for row in rows:
        if row["catalyst"] != catalyst or row.get("facet", "") != facet:
            continue
        status = row.get("converged", "true").lower()
        if status not in ("true", "false", "1", "0", "yes", "no", ""):
            raise ValueError("Invalid converged value: {}".format(status))
        if status in ("false", "0", "no"):
            skipped += 1
            continue
        species = "+".join(aliases.get(p, p) for p in row["species"].split("+"))
        energy = finite(row["free_energy"], "free_energy of " + species)
        entry = dict(row, free_energy=energy)
        if species not in selected or (energy, row.get("structure_id", "")) < (
                selected[species]["free_energy"], selected[species].get("structure_id", "")):
            selected[species] = entry
    if not selected:
        raise ValueError("No converged energies for catalyst={!r}, facet={!r}".format(catalyst, facet))
    return selected, skipped


def import_legacy_csv(source, destination, column="free_energy", facet=""):
    """Import old headerless 'Ga_CO,-1.23' without guessing barrier semantics."""
    import csv
    if column not in ("free_energy", "energy"):
        raise ValueError("column must be free_energy or energy")
    rows, seen = [], set()
    with open(source, encoding="utf-8-sig", newline="") as stream:
        for line, values in enumerate(csv.reader(stream), 1):
            if not values:
                continue
            if len(values) != 2 or "_" not in values[0]:
                raise ValueError("Line {}: expected catalyst_species,value".format(line))
            catalyst, species = values[0].strip().split("_", 1)
            if not catalyst or not species or (catalyst, species) in seen:
                raise ValueError("Empty or duplicate key on line {}".format(line))
            seen.add((catalyst, species))
            rows.append({"catalyst": catalyst, "facet": facet, "species": species,
                         column: finite(values[1], "line {}".format(line))})
    write_csv(destination, rows, ["catalyst", "facet", "species", column])
    return len(rows)

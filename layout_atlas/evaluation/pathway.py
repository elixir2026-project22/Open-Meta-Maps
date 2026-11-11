"""Pathway graph from SBML (species + reactions as nodes).

Adapted from layout-benchmark `graph.py` (MIT). Kept lean: only what evaluation needs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import libsbml
import networkx as nx

SPECIES = "species"
REACTION = "reaction"
REACTANT, PRODUCT, MODIFIER = "reactant", "product", "modifier"
SPECIES_SIZE = (60.0, 24.0)
REACTION_SIZE = (10.0, 10.0)


@dataclass
class Pathway:
    name: str
    graph: nx.DiGraph
    compartments: dict[str, str] = field(default_factory=dict)
    source: Path | None = None

    @property
    def species(self) -> list[str]:
        return [n for n, d in self.graph.nodes(data=True) if d["kind"] == SPECIES]

    @property
    def reactions(self) -> list[str]:
        return [n for n, d in self.graph.nodes(data=True) if d["kind"] == REACTION]

    def undirected(self) -> nx.Graph:
        g = nx.Graph()
        g.add_nodes_from(self.graph.nodes(data=True))
        g.add_edges_from((u, v) for u, v in self.graph.edges() if u != v)
        return g


def _add_species(g, sid, label, compartment):
    w, h = SPECIES_SIZE
    g.add_node(sid, kind=SPECIES, label=label or sid, compartment=compartment, width=w, height=h)


def _add_reaction(g, rid, label, compartment=None):
    w, h = REACTION_SIZE
    g.add_node(rid, kind=REACTION, label=label or rid, compartment=compartment, width=w, height=h)


def read_sbml(path: str | Path, name: str | None = None) -> Pathway:
    path = Path(path)
    doc = libsbml.readSBMLFromFile(str(path))
    if doc.getNumErrors(libsbml.LIBSBML_SEV_FATAL) or doc.getModel() is None:
        raise ValueError(f"{path}: not a readable SBML model")
    model = doc.getModel()
    g = nx.DiGraph()
    compartments = {c.getId(): c.getName() or c.getId() for c in model.getListOfCompartments()}
    for s in model.getListOfSpecies():
        _add_species(g, s.getId(), s.getName(), s.getCompartment())
    for r in model.getListOfReactions():
        rid = r.getId()
        _add_reaction(g, rid, r.getName(), r.getCompartment() if r.isSetCompartment() else None)
        for sr in r.getListOfReactants():
            g.add_edge(sr.getSpecies(), rid, role=REACTANT)
        for sr in r.getListOfProducts():
            g.add_edge(rid, sr.getSpecies(), role=PRODUCT)
        for sr in r.getListOfModifiers():
            if not g.has_edge(sr.getSpecies(), rid):
                g.add_edge(sr.getSpecies(), rid, role=MODIFIER)
    g.remove_nodes_from([n for n in list(g) if g.degree(n) == 0])
    return Pathway(name=name or path.stem, graph=g, compartments=compartments, source=path)

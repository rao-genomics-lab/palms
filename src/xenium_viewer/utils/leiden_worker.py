"""Subprocess worker for Leiden clustering.

Zero top-level imports — this module must be importable quickly in a spawned
subprocess (avoids pulling in napari/scanpy/Qt in the child process).
"""


def run_leiden(data, indices, indptr, n, resolution, seed=42):
    """Run leidenalg in a subprocess. All heavy imports are local."""
    import scipy.sparse as sp
    import igraph as ig
    import leidenalg

    conn = sp.csr_matrix((data, indices, indptr), shape=(n, n))
    cx = conn.tocoo()
    g = ig.Graph(
        n=n,
        edges=list(zip(cx.row.tolist(), cx.col.tolist())),
        directed=True,
    )
    g.es['weight'] = cx.data.tolist()
    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        weights='weight',
        resolution_parameter=resolution,
        seed=seed,
    )
    return partition.membership  # list[int]

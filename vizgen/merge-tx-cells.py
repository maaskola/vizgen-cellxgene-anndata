#!/usr/bin/env python

import os
import os.path as osp
from typing import Optional
import click
import h5py
import numpy as np
import pandas as pd
# import dask.dataframe as dd
import cudf
import cupy as cp
import cuspatial

import matplotlib.pyplot as plt

def pd_read_tx(path, nrows=None):
    dtype = {
        'barcode_id': np.uintc,
        'global_x': np.float64,
        'global_y': np.float64,
        'global_z': np.uintc,
        # 'global_z': np.float64,
        'x': np.float64,
        'y': np.float64,
        'fov': np.uintc,
        'gene': 'category',
        'transcript_id': 'category',
    }
    return pd.read_csv(path, dtype=dtype, index_col=0, header=0, nrows=nrows)

def pd_read_cell_metadata(path):
    dtype = {
        'fov': np.uintc,
        'volume': np.float64,
        'center_x': np.float64,
        'center_y': np.float64,
        'min_x': np.float64,
        'max_x': np.float64,
        'min_y': np.float64,
        'max_y': np.float64,
    }
    return pd.read_csv(path, dtype=dtype, index_col=0, header=0)

def read_cell_boundaries_file(path, z_level):
    component = 'zIndex_{}/p_0/coordinates'.format(int(z_level))
    # print(component)
    h5py_file = h5py.File(path)
    return {
        cell: np.squeeze(np.array(content[component]))
        for cell, content in h5py_file['featuredata'].items()
    }

def determine_extents(xs, ys, max_depth):
    x_min = xs.min()
    x_max = xs.max()
    y_min = ys.min()
    y_max = ys.max()
    scale = max(x_max - x_min, y_max - y_min) // (1 << max_depth)
    return x_min, x_max, y_min, y_max, scale

def extract_polygons(polygons):
    n = len(polygons)
    polys = polygons.values()
    poly_offsets = cudf.Series(cp.arange(n), index=polygons.keys())
    ring_offsets = cp.cumsum(cp.array([0] + [len(poly) for poly in polys][:-1]))
    poly_points_x = cp.array(np.concatenate([poly[:, 0] for poly in polys]))
    poly_points_y = cp.array(np.concatenate([poly[:, 1] for poly in polys]))
    return poly_offsets, ring_offsets, poly_points_x, poly_points_y

def visualize(
        polygons,
        x,
        y,
        z_level: int,
        dir_path: Optional[str],
        fov: int,
        point_scale: float = 2.0,
):
    plt.figure(figsize=(32, 32))
    plt.axis('equal')
    for cell_id in polygons:
        coords = polygons[cell_id].squeeze()
        plt.fill(coords[:, 0], coords[:, 1])
    plt.scatter(x, y, s=point_scale, alpha=1.0, c='black')
    # plt.show()
    if osp.exists(dir_path):
        if not osp.isdir(dir_path):
            raise RuntimeError('Visualization output path exists and is not a directory: ' + dir_path)
    else:
        os.mkdir(dir_path)

    plt.savefig(osp.join(dir_path, 'fov-{}-z{}.png'.format(fov, z_level)))
    plt.close()

def intersection(polygons, xs, ys, max_depth=3, min_size=50):
    x_min, x_max, y_min, y_max, scale = determine_extents(xs, ys, max_depth)
    poly_offsets, ring_offsets, poly_points_x, poly_points_y = extract_polygons(polygons)
    # print([x_min, x_max, y_min, y_max, scale])
    keys_to_points, quadtree = cuspatial.quadtree_on_points(
        xs,
        ys,
        x_min,
        x_max,
        y_min,
        y_max,
        scale,
        max_depth,
        min_size,
    )
    poly_bounding_boxes = cuspatial.polygon_bounding_boxes(
        poly_offsets,
        ring_offsets,
        poly_points_x,
        poly_points_y,
    )
    poly_quad_pairs = cuspatial.join_quadtree_and_bounding_boxes(
        quadtree,
        poly_bounding_boxes,
        x_min,
        x_max,
        y_min,
        y_max,
        scale,
        max_depth,
    )
    res_df = cuspatial.quadtree_point_in_polygon(
        poly_quad_pairs,
        quadtree,
        keys_to_points,
        xs,
        ys,
        poly_offsets,
        ring_offsets,
        poly_points_x,
        poly_points_y,
    )
    # assert each point intersects with at most one polygon
    uniq = res_df['point_index'].unique()
    n_uniq = len(uniq)
    assert(n_uniq == res_df.shape[0])
    idxs = cp.full((len(xs),), fill_value=-1)
    idxs[res_df['point_index']] = res_df['polygon_index']
    return cp.asnumpy(idxs), n_uniq

def intersect_tx_cells(
        fov: int,
        z_level: int,
        tx: pd.DataFrame,
        boundaries_dir: str,
        viz_dir: Optional[str],
        *args,
        **kw_args
):
    # print('doing FOV {} with {} entries'.format(fov, len(tx)))
    polygons = read_cell_boundaries_file(
        osp.join(
            boundaries_dir,
            'feature_data_{}.hdf5'.format(fov)
        ),
        z_level,
    )
    tx_x = tx['global_x']
    tx_y = tx['global_y']
    if viz_dir is not None:
        visualize(polygons, tx_x, tx_y, z_level, dir_path=viz_dir, fov=fov)
    intersected, n_uniq = intersection(polygons, tx_x, tx_y, *args, **kw_args)
    print('FOV {}: {} of {} tx molecules in cells: {} %'.format(
        fov,
        n_uniq,
        tx.shape[0],
        n_uniq / tx.shape[0] * 100,
    ))
    return intersected

def get_angle(direction):
    # TODO
    pass

def do_angles(tx, metadata):
    tx_coords = tx[['x', 'y']]
    center_coords = metadata.iloc(tx['cell'])[['center_x', 'center_y']]
    return get_angle(tx_coords - center_coords)

@click.command()
@click.argument('tx_path', type=click.Path(exists=True, dir_okay=False))
@click.argument('cell_metadata_path', type=click.Path(exists=True, dir_okay=False))
@click.argument('cell_boundaries_dir_path', type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.argument('out_path', type=click.Path(exists=False))
@click.option('--nrows', type=int, help='number of transcripts to consider')
@click.option('--max-depth', default=3, type=int, show_default=True, help='maximum quadtree depth')
@click.option('--min-size', default=50, type=int, show_default=True, help='minimum number of points for a non-leaf quadtree node')
@click.option('--viz-dir', type=click.Path(dir_okay=True, file_okay=False), show_default=True, help='directory where FOVs are visualized')
def run(
        tx_path: str,
        cell_metadata_path: str,
        cell_boundaries_dir_path: str,
        out_path: str,
        nrows: Optional[int],
        max_depth: int,
        min_size: int,
        viz_dir: str,

):
    metadata_grouped = pd_read_cell_metadata(cell_metadata_path).groupby('fov')
    tx_df = pd_read_tx(tx_path, nrows=nrows)
    # df = dd.read_csv(tx_path)
    blocks = []
    for fov, fov_group in tx_df.groupby('fov'):
        print('FOV {}: before filtering blank genes: {}'.format(fov, len(fov_group)))
        fov_group = fov_group[~fov_group['gene'].str.startswith('Blank')]
        print('FOV {}: after filtering blank genes: {}'.format(fov, len(fov_group)))
        for z_level, z_group in fov_group.groupby('global_z'):
            # print('FOV {}: z-level {}: {}'.format(fov, z_level, len(z_group)))
            z_group['cell'] = intersect_tx_cells(
            # intersect_tx_cells(
                fov=fov,
                z_level=z_level,
                tx=z_group,
                boundaries_dir=cell_boundaries_dir_path,
                viz_dir=viz_dir,
                max_depth=max_depth,
                min_size=min_size,
            )
            # breakpoint()
            # z_group['angle'] = do_angles(z_group, metadata_grouped[fov])
            blocks += [z_group]
        # breakpoint()
    joined = pd.concat(blocks)
    joined.to_csv(out_path)

if __name__ == '__main__':
    run()

import anndata as ad
import pandas as pd
import click

@click.command()
@click.argument('csv-path')
@click.option('--out-path', default='cell_by_gene.h5ad', show_default=True)
@click.option('--metadata-path', default='cell_metadata.csv', show_default=True)
def run(csv_path: str,
        out_path: str,
        metadata_path: str,
        ):
    adata = ad.read_csv(csv_path, first_column_names=True)
    print('loaded {}'.format(adata))
    if metadata_path:
        metadata = pd.read_csv(metadata_path, header=0, index_col=0)
        assert adata.n_obs == len(metadata)
        print('loaded metadata {}'.format(metadata_path))
        for column in metadata.columns:
            adata.obs[column] = metadata[column]
    print('writing to {}'.format(out_path))
    adata.write(out_path, compression="gzip")

if __name__ == '__main__':
    run()

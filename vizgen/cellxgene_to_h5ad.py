import anndata as ad
import click

@click.command()
@click.argument('csv-path')
@click.option('--out-path', default='cell_by_gene.h5ad', show_default=True)
def run(csv_path: str,
        out_path: str,
        ):
    adata = ad.read_csv(csv_path, first_column_names=True)
    print('loaded {}'.format(adata))
    print('writing to {}'.format(out_path))
    adata.write(out_path, compression="gzip")

if __name__ == '__main__':
    run()

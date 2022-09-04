import os
import json
import logging
import logging.config as log_config
import click
import builder
from builder.databases.parsers import parsers
from joblib import Parallel, delayed

config_dir = os.path.dirname(os.path.abspath(builder.__file__))

# print(sys.path)
with open(os.path.join(config_dir, "log.config")) as f:
    config = json.load(f)
    log_config.dictConfig(config)

logger = logging.getLogger(__name__)


@click.group()
def database():
    pass


def parse_database(import_directory, database_directory, database, download=True, skip=True):
    stats = set()
    Parser = parsers.get(database, None)
    if Parser:
        parser = Parser(import_directory, database_directory, download, skip)
        stats = parser.build_stats()
    return stats


@database.command()
@click.option('--db-dir', '-d', required=True,
              type=click.Path(exists=True, dir_okay=True),
              help="The directory which saved the downloaded database files.")
@click.option('--output-dir', '-o', required=True,
              type=click.Path(exists=True, dir_okay=True),
              help="The directory which saved the graph files.")
@click.option('--database', required=True, type=click.Choice(parsers.keys()),
              help="Which database?", multiple=True)
@click.option('--n-jobs', '-n', required=False,
              help="Hom many jobs?", default=4)
@click.option('--download/--no-download', default=False)
@click.option('--skip/--no-skip', default=True)
def make_database(output_dir, db_dir, database, download, n_jobs, skip):
    all_databases = database
    valid_databases = list(
        filter(lambda database: database in parsers.keys(), all_databases))
    invalid_databases = list(
        filter(lambda database: database not in parsers.keys(), all_databases))
    if len(invalid_databases) > 0:
        logger.warn("%s databases (%s) is not valid, skip them.",
                    len(invalid_databases), invalid_databases)
    stats = Parallel(n_jobs=n_jobs)(delayed(parse_database)(
        output_dir, db_dir, database, download, skip) for database in valid_databases)
    allstats = {val if type(
        sublist) == set else sublist for sublist in stats for val in sublist}
    return allstats


if __name__ == "__main__":
    main = click.CommandCollection(sources=[database])
    main()

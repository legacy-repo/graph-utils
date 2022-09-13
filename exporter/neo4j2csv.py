#!/usr/bin/env python3

import re
import os
import csv
import click
import json
from py2neo import Graph


class CustomRecord:
    def __init__(self, record, identity, show_attrs) -> None:
        self.record = record
        self.show_attrs = show_attrs
        self.identity = identity
        self.rel_record = record.get('r')
        if self.rel_record is None:
            raise Exception("You need to return r in your cypher clause.")

        self.node_record = record.get('n')
        if self.node_record is None:
            raise Exception("You need to return n in your cypher clause.")

    def format(self, inst):
        if type(inst) == str:
            return inst
        elif type(inst) == list:
            return ';'.join(inst)
        elif type(inst) == dict:
            return json.dump(inst)

    def format_items(self, items):
        output = {}
        for i in items:
            output[i[0]] = self.format(i[1])
        return output

    def get_relationship(self):
        raise NotImplemented

    def get_node(self):
        raise NotImplemented


class GephiRecord(CustomRecord):
    def __init__(self, record, identity, show_attrs) -> None:
        super().__init__(record, identity, show_attrs)

    def get_relationship(self):
        output = {
            "Source": self.format(self.rel_record.start_node.get(self.identity, None)),
            "Target": self.format(self.rel_record.start_node.get(self.identity, None)),
            "Type": self.format(type(self.rel_record).__name__)
        }
        return output

    def get_node(self):
        output = {
            "Id": self.format(self.node_record.get(self.identity, None)),
            "Label": self.format(self.node_record.labels),
        }
        if self.show_attrs:
            output.update(**self.format_items(self.node_record.items()))
        return output


class DglkeRecord(CustomRecord):
    def __init__(self, record, identity, show_attrs) -> None:
        super().__init__(record, identity, show_attrs)

    def get_relationship(self):
        output = {
            "start": self.format(self.rel_record.start_node.get(self.identity, None)),
            "end": self.format(self.rel_record.start_node.get(self.identity, None)),
            "type": self.format(type(self.rel_record).__name__)
        }
        return output

    def get_node(self):
        output = {
            "id": self.format(self.node_record.get(self.identity, None)),
            "label": self.format(self.node_record.labels),
        }
        if self.show_attrs:
            output.update(**self.format_items(self.node_record.items()))
        return output


def format_record(record, identity, show_attrs=True, output_format="gephi"):
    if output_format == 'gephi':
        gephi_record = GephiRecord(record, identity, show_attrs)
        node = gephi_record.get_node()
        relationship = gephi_record.get_relationship()
    else:
        dglke_record = DglkeRecord(record, identity, show_attrs)
        node = dglke_record.get_node()
        relationship = dglke_record.get_relationship()

    return node, relationship


def write2csv(fhandler, odict: dict, is_first=False):
    w = csv.DictWriter(fhandler, odict.keys())
    if is_first:
        w.writeheader()
    w.writerow(odict)


def write_batch(fname, dicts):
    with open(fname, 'w+') as f:
        for idx, i in enumerate(dicts):
            write2csv(f, i, True if idx == 0 else False)


@click.command()
@click.option('--uri', required=True, help='Where neo4j instance(such as `g.example.com:7687)?`')
@click.option('--database', '-d', required=True, help='Which database?')
@click.option('--username', default='neo4j', help='The username of neo4j instance.')
# @click.option('--password', prompt=True, hide_input=True, confirmation_prompt=True, help='The password of neo4j instance.')
@click.password_option(help='The password of neo4j instance.')
@click.option('--query-str', help='Cypher query string. such as `MATCH (n) OPTIONAL MATCH (n)-[r]-() RETURN n, r SKIP %s LIMIT %s;`')
@click.option('--identity', default='id', help='What name is used to identify a node?')
@click.option('--show-attrs/--no-show-attrs', default=True, help='Whether output all the attributes?')
@click.option('--output-dir', '-o', default='.', help='The output directory')
@click.option('--output-format', default='gephi', type=click.Choice(['gephi', 'dglke']), help='The format of output file?')
@click.option('--batch', '-b', default=0, help='which batch do you want to start?')
@click.option('--batch-size', '-s', default=10, type=click.IntRange(0, 1000, clamp=True), help='How many records do you want to output in one batch?')
def export(uri, database, username, password, query_str, show_attrs, output_dir, identity, output_format, batch, batch_size):
    """Export subgraph to csv."""
    print("Connect to %s" % uri)
    print("Output files are located in %s" % os.path.abspath(output_dir))
    if username and password:
        graph = Graph("neo4j://%s" % uri, name=database,
                      auth=(username, password))
    else:
        graph = Graph("bolt+s://%s" % uri, name=database)

    if not re.match(r'.*RETURN n,r', query_str):
        raise Exception("Your cyper clause need to return nodes and relationships by using RETURN, such as `RETURN n,r`")

    if query_str:
        query = graph.run("%s SKIP %s LIMIT %s;" % (query_str, batch * batch_size, batch_size))
        obatch = [format_record(i, identity, show_attrs, output_format) for i in query]
    else:
        query = graph.run("MATCH (n) OPTIONAL MATCH (n)-[r]-() RETURN n, r SKIP %s LIMIT %s;" % (batch * batch_size, batch_size))
        obatch = [format_record(i, identity, show_attrs, output_format) for i in query]

    write_batch(os.path.join(output_dir, "nodes.csv"), list(map(lambda x: x[0], obatch)))
    write_batch(os.path.join(output_dir, "relationship.csv"), list(map(lambda x: x[1], obatch)))
        


if __name__ == '__main__':
    export()

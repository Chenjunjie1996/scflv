#!/usr/bin/env python

import collections
import csv
import gzip
import os
import re
import sys

PATTERN = re.compile(r'(\S+?)\s*"(.*?)"')
gtf_row = collections.namedtuple( 'gtf_row', 'seqname source feature start end score strand frame attributes' )

def generic_open(file_name, *args, **kwargs):
    if file_name.endswith('.gz'):
        file_obj = gzip.open(file_name, *args, **kwargs)
    else:
        file_obj = open(file_name, *args, **kwargs)
    return file_obj


class GtfParser:
    def __init__(self, gtf_fn):
        self.gtf_fn = gtf_fn
        self.gene_id = []
        self.gene_name = []
        self.id_name = {}
        self.id_strand = {}

    def get_properties_dict(self, properties_str):
        """
        allow no space after semicolon
        """

        if isinstance(properties_str, dict):
            return properties_str

        properties = collections.OrderedDict()
        attrs = properties_str.split(';')
        for attr in attrs:
            if attr:
                m = re.search(PATTERN, attr)
                if m:
                    key = m.group(1).strip()
                    value = m.group(2).strip()
                    properties[key] = value

        return properties

    def gtf_reader_iter(self):
        """
        Yield:
            row: list
            gtf_row
        """
        with generic_open(self.gtf_fn, mode='rt') as f:
            reader = csv.reader(f, delimiter='\t')
            for i, row in enumerate(reader, start=1):
                if len(row) == 0:
                    continue
                if row[0].startswith('#'):
                    yield row, None
                    continue

                if len(row) != 9:
                    sys.exit(f"Invalid number of columns in GTF line {i}: {row}\n")

                if row[6] not in ['+', '-']:
                    sys.exit(f"Invalid strand in GTF line {i}: {row}\n")

                seqname = row[0]
                source  = row[1]
                feature = row[2]
                # gff/gtf is 1-based, end-inclusive
                start   = int(row[3])
                end     = int(row[4])
                score   = row[5]
                strand  = row[6]
                frame  = row[7]
                attributes = self.get_properties_dict(row[8])

                yield row, gtf_row( seqname, source, feature, \
                    start, end, score, strand, frame, \
                    attributes )

def filter_gtf(gtf_fn, out_fn, allow):
    """
    Filter attributes

    Args:
        allow: {
            "gene_biotype": set("protein_coding", "lncRNA")
        }
    """
    sys.stderr.write("Writing GTF file...\n")
    gp = GtfParser(gtf_fn)
    n_filter = 0

    with open(out_fn, 'w') as f:
        # quotechar='' is not allowed since python3.11
        writer = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar=None)
        for row, grow in gp.gtf_reader_iter():
            if not grow:
                writer.writerow(row)
                continue


            remove = False
            if allow:
                for key, value in grow.attributes.items():
                    if key in allow and value not in allow[key]:
                        remove = True
                        break

            if not remove:
                writer.writerow(row)
            else:
                n_filter += 1
    return n_filter


if __name__ == '__main__':
    # args: gtf, attributes
    gtf_fn = sys.argv[1]
    attributes = sys.argv[2]
    out_fn = os.path.basename(gtf_fn).replace('.gtf', '.filtered.gtf')

    allow = {}
    for attr_str in attributes.split(';'):
        if attr_str:
            attr, val = attr_str.split('=')
            val = set(val.split(','))
            allow[attr] = val

    n_filter = filter_gtf(gtf_fn, out_fn, allow)
    sys.stdout.write(f"Filtered {n_filter} lines\n")
    log_file = "gtf_filter.log"
    with open(log_file, 'w') as f:
        f.write(f"Filtered lines: {n_filter}\n")
        f.write(f"Attributes: {attributes}\n")
        f.write(f"Output file: {out_fn}\n")



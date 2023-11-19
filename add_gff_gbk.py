import argparse
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF

def arg_parser():
    parser = argparse.ArgumentParser(description="Add GFF information to GenBank file")
    parser.add_argument("-g", "--gff", help="GFF file", required=True)
    parser.add_argument("-f", "--gbk", help="GenBank file", required=True)
    parser.add_argument("-o", "--output", help="Output file", required=True)
    return parser.parse_args()

def add_gff_to_gbk(gff_file, gbk_file, output_file):
    # Read the GenBank file
    gbk_records = SeqIO.to_dict(SeqIO.parse(gbk_file, "genbank"))

    # Initialize a dictionary to store GFF information
    gff_data = {}

    # Parse the GFF file and store the relevant information
    with open(gff_file, "r") as gff_handle:
        for rec in GFF.parse(gff_handle):
            for feature in rec.features:
                if feature.type == "rRNA":
                    locus_tag = feature.qualifiers["Name"][0]
                    product = feature.qualifiers["product"][0]
                    start = feature.location.start.position
                    end = feature.location.end.position
                    strand = "+" if feature.location.strand == 1 else "-"

                    if rec.id not in gff_data:
                        gff_data[rec.id] = []

                    gff_data[rec.id].append({
                        "locus_tag": locus_tag,
                        "product": product,
                        "start": start,
                        "end": end,
                        "strand": strand
                    })

    # Update the GenBank records with GFF information
    for record_id, gff_features in gff_data.items():
        if record_id in gbk_records:
            for feature_info in gff_features:
                gbk_record = gbk_records[record_id]
                new_feature = SeqFeature()  # Use SeqFeature from Bio.SeqFeature
                new_feature.location = FeatureLocation(  # Use FeatureLocation from Bio.SeqFeature
                    start=feature_info["start"],
                    end=feature_info["end"],
                    strand=1 if feature_info["strand"] == "+" else -1
                )
                new_feature.type = "rRNA"
                new_feature.qualifiers["locus_tag"] = feature_info["locus_tag"]
                new_feature.qualifiers["product"] = feature_info["product"]

                gbk_record.features.append(new_feature)

    # Write the updated GenBank records to a new file
    with open(output_file, "w") as output_handle:
        for record_id, gbk_record in gbk_records.items():
            SeqIO.write(gbk_record, output_handle, "genbank")
if __name__ == "__main__":
    args = arg_parser()
    add_gff_to_gbk(args.gff, args.gbk, args.output)

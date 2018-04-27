# guide for downloading allen brain atlas SVG images
# https://gist.github.com/lydiang/12d04cd4ed750f4b24cfc72d67d14387

import argparse
import os
import requests
from pathlib import Path

output_directory = Path('./data/svg/slices')

# Specify downsample factor
downsample = 6

# Parameters from: http://help.brain-map.org/display/api/Atlas+Drawings+and+Ontologies
# "Atlas ID"
adult_human = 265297125 # Human Brain Atlas Guide
fetal15 = 138322603 # 15pcw fetal human brain
fetal21 = 3 # 21pcw fetal human brain
fetal21_brainstem = 287730656 # 21pcw fetal human brainstem

# "GraphicGroupLabels"
fetal21_labels = [31,113753816,141667008] #fetal21
fetal15_labels = [31,141667008] #fetal15,  fetal21_brainstem
adult_labels = [265297119, 266932194, 266932196, 266932197] #adult

atlas_info = {'atlas_id': {'adult': adult_human,
                           'fetal15': fetal15,
                           'fetal21': fetal21,
                           'fetal21_brainstem': fetal21_brainstem},
              'group_labels': {'adult': adult_labels,
                               'fetal15': fetal15_labels,
                               'fetal21': fetal21_labels,
                               'fetal21_brainstem': fetal15_labels}}

def get_image_info(atlas_id):

    # RMA query to find images for atlas
    query_url = "http://api.brain-map.org/api/v2/data/query.json?criteria=model::AtlasImage"
    query_url += ",rma::criteria,[annotated$eqtrue]"
    query_url += ",atlas_data_set(atlases[id$eq%d])" % (atlas_id)
    query_url += ",rma::options[order$eq'sub_images.section_number'][num_rows$eqall]"

    r = requests.get(query_url)
    image_info = r.json()['msg']
    return image_info


def write_images(image_info, downsample):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # loop through each image
    for image in image_info:

        print(image['section_number'])
        # download svg
        svg_url = "http://api.brain-map.org/api/v2/svg_download/%d?" % (image['id'])
        svg_url += "groups=%s" % (",".join([str(g) for g in group_labels]))
        svg_url += "&downsample=%d" % (downsample)
        svg_path = os.path.join(output_directory, '%04d_%d.svg' % (image['section_number'], image['id']))
        r = requests.get(svg_url)
        with open(svg_path, 'wb') as f:
            f.write(r.content)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("atlas", type=str,
                        help="Select atlas: adult, fetal15, fetal21, fetal21_brainstem")
    parser.add_argument("--downsample", '-d', type=int, default=6,
                        help="define downsampling factor (default is 6)")
    args = parser.parse_args()

    atlas_id = atlas_info['atlas_id'][args.atlas]
    group_labels = atlas_info['group_labels'][args.atlas]

    image_info = get_image_info(atlas_id)
    output_directory = output_directory / args.atlas
    print('There are {} images. Saving to disk at {}'.format(len(image_info), output_directory))
    write_images(image_info, args.downsample)

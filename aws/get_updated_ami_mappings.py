#!/usr/bin/env python

import sys
import os
import argparse
import boto3


def _parse_arguments(desc, theargs):
    """Parses command line arguments using argparse
    """
    help_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_formatter)
    parser.add_argument("--ownerid", default='679593333241',
                        help="Owner id to pass to search " +
                             "(default 679593333241)")
    parser.add_argument('--namefilter',
                        default='*1805_01*',
                        help='Find only AMI image with this string in name' +
                             ' (default: *1805_01*)')
    parser.add_argument("--productcode", 
                        default='aw0evgkw8e5c1q413zgy5pjce',
                        help="product code to filter with" +
                             "(default aw0evgkw8e5c1q413zgy5pjce)")
    parser.add_argument('--profile',
                        default=None,
                        help='AWS profile to load from credentials. default none')
    return parser.parse_args(theargs)

def _get_ami_mapping(theargs):
    """Returns a string containing ami mapping
    """
    import pprint
    pp = pprint.PrettyPrinter(indent=2)
    mapstr = ''
    if theargs.profile is not None:
        boto3.setup_default_session(profile_name=theargs.profile)

    ec2 = boto3.client('ec2')
    response = ec2.describe_regions()
    for region in response['Regions']:
        rname = region['RegionName']
        sys.stdout.write('Running query in region: ' + rname + '\n')
        ec2 = boto3.client('ec2', region_name=rname)
        resp = ec2.describe_images(Owners=[theargs.ownerid],
                                   Filters=[{'Name': 'product-code',
                                             'Values': [theargs.productcode]},
                                             {'Name': 'virtualization-type',
                                             'Values': ['hvm']},
                                             {'Name': 'name',
                                             'Values': [theargs.namefilter]}])
        pp.pprint(resp)
        for image in resp['Images']:
           
            mapstr += ('           "' + rname + '"   : {"AMI" : "' +
                       image['ImageId'] + '"},\n')
    sys.stdout.write('\n\n Below is json fragment that can ' + 
                     'go in "RegionMap" of cloud formation template\n\n')
    return mapstr


def main(arglist):
    desc = """
              This script uses AWS boto library to query for AMI
              images that match the owner and name filter
              passed in via --ownerid and --namefilter flags for
              this tool.  The output is json fragment that can
              be put in the "Mappings" => "RegionMap" section

           """
    theargs = _parse_arguments(desc, sys.argv[1:])
    sys.stdout.write('Querying AWS: \n')
    sys.stdout.write(_get_ami_mapping(theargs))


if __name__ == '__main__': # pragma: no cover
    sys.exit(main(sys.argv))

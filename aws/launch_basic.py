#!/usr/bin/env python

import sys
import os
import argparse
import boto3
from ipify import get_ip


def _parse_arguments(desc, theargs):
    """Parses command line arguments using argparse
    """
    help_formatter = argparse.RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=help_formatter)
    parser.add_argument('--template', required=True,
                        help='CloudFormation template file to use')
    parser.add_argument('--region', default='us-east-2',
                        help="Region to use" +
                             "(default us-east-2)")
    parser.add_argument('--name', default='USERNAMEstack',
                        help='Stack name to use')
    parser.add_argument('--profile', default=None,
                        help='AWS profile to load from credentials. default none')
    parser.add_argument('--keypairname', default='id_rsa',
                        help='AWS EC2 KeyPair Name')
    parser.add_argument('--instancetype', default='c5.xlarge',
                        help='Instance type to launch (default c5.xlarge')
    parser.add_argument('--disksize', default='1000',
                        help='Disk Size in gigabytes (default 1000)')
    parser.add_argument('--accesslocation', default='',
                        help='ip4 CIDR to denote ip address(s) to allow '
                             'http/ssh access to EC2 instance. (default is ip '
                             'address of machine running this script')
    return parser.parse_args(theargs)


def _launch_cloudformation(theargs):
    """Launches cloud formation
    """
    if theargs.profile is not None:
        boto3.setup_default_session(profile_name=theargs.profile)

    cloudform = boto3.client('cloudformation', region_name=theargs.region)
    template = theargs.template
    with open(template, 'r') as f:
        template_data = f.read()

    if theargs.accesslocation is None or theargs.accesslocation is '':
        theargs.accesslocation = str(get_ip()) + '/32'

    params = [
        {
            'ParameterKey': 'KeyName',
            'ParameterValue': theargs.keypairname
        },
        {
            'ParameterKey': 'InstanceType',
            'ParameterValue': theargs.instancetype
        },
        {
            'ParameterKey': 'DiskSize',
            'ParameterValue': theargs.disksize
        },
        {
            'ParameterKey': 'AccessLocation',
            'ParameterValue': theargs.accesslocation
        }
    ]

    tags = [
        {
            'Key': 'Name',
            'Value': theargs.name
        }
    ]

    resp = cloudform.create_stack(
        StackName=theargs.name,
        TemplateBody=template_data,
        Parameters=params,
        TimeoutInMinutes=25,
        Tags=tags
    )
    return str(resp)


def main(arglist):
    desc = """
              Launches CloudFormation template
           """
    theargs = _parse_arguments(desc, sys.argv[1:])
    sys.stdout.write('Contacting AWS: \n')
    sys.stdout.write(_launch_cloudformation(theargs))


if __name__ == '__main__': # pragma: no cover
    sys.exit(main(sys.argv))

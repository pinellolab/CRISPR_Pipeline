#!/usr/bin/env python

import argparse


def apply_spacer_logic(raw_technology: str) -> str:
    parts = raw_technology.split(':')
    if len(parts) < 3:
        return raw_technology

    transcript = parts[2].split(',')
    if len(transcript) < 3:
        return raw_technology

    transcript[0] = transcript[0] or '1'
    transcript[1] = '0'
    transcript[2] = '0'
    parts[2] = ','.join(transcript)
    return ':'.join(parts)


def main():
    parser = argparse.ArgumentParser(description='Prepare guide chemistry and workflow settings for kb count.')
    parser.add_argument('--technology-file', required=True, help='Path to a file containing the raw chemistry string.')
    parser.add_argument('--is-10x3v3', required=True, help='Whether to force 10XV3/kite:10xFB behavior.')
    parser.add_argument('--spacer-tag', default='', help='Optional spacer tag prepended to guides.')
    args = parser.parse_args()

    with open(args.technology_file) as handle:
        raw_technology = handle.read().strip()

    is_10x3v3 = str(args.is_10x3v3).lower() == 'true'
    spacer_tag = args.spacer_tag or ''

    if is_10x3v3:
        technology = '10XV3'
        workflow = 'kite:10xFB'
    else:
        technology = apply_spacer_logic(raw_technology) if spacer_tag else raw_technology
        workflow = 'kite'

    with open('kb_technology.txt', 'w') as handle:
        handle.write(f'{technology}\n')

    with open('kb_workflow.txt', 'w') as handle:
        handle.write(f'{workflow}\n')


if __name__ == '__main__':
    main()
import argparse
import re
import sys

parser = argparse.ArgumentParser(
    prog='indent_comments',
    description='Ensures that Doxygen comments are indented consistently')

parser.add_argument('filename', nargs="+")
parser.add_argument('-f', '--fix', action='store_true')
parser.add_argument('-i', '--inplace', action='store_true')

args = parser.parse_args()
if args.inplace and not args.fix:
    print("-i requires -f", file=sys.stderr)
    sys.exit(2)

broken = False

for file in args.filename:
    with open(file, "a+t" if args.inplace else "rt") as fh:
        indent = None
        lines = fh
        out = sys.stdout
        if args.inplace:
            fh.seek(0)
            lines = fh.readlines()
            fh.seek(0)
            fh.truncate()
            out = fh
        for nr, line in enumerate(lines, start=1):
            if indent is not None:
                m = re.match("^([ \t]*)(\\*.*$)", line)
                if not m:
                    print(
                        f"{file}:{nr} {line!r} should be inside multi-line doxygen comment, but doesn't start with '*'",
                        file=sys.stderr)
                    broken = True
                    if args.fix:
                        print(line, end="", file=out)
                elif args.fix:
                    print(indent, m.group(2), file=out)  # inserts single space between args
                elif indent + " " != m.group(1):
                    print(f"{file}:{nr} has indent {m.group(1)!r} but should have {indent + ' '!r}")
                    broken = True
                if "*/" in line:
                    indent = None
            else:
                if args.fix:
                    print(line, end="", file=out)
                m = re.match("^([ \t]*)/\\*\\*", line)
                if not m or "*/" in line: continue
                indent = m.group(1)

if broken:
    sys.exit(1)

import argparse
import fileinput as fi
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
indent = None

for line in fi.input(files=args.filename, inplace=args.inplace):
    if fi.isfirstline():
        indent = None
    if indent is not None:
        m = re.match("^([ \t]*)(\\*.*$)", line)
        if not m:
            print(
                f"{fi.filename()}:{fi.filelineno()} {line!r} should be inside multi-line doxygen comment, but doesn't start with '*'",
                file=sys.stderr)
            broken = True
            if args.fix:
                print(indent, "*", line.strip())

        elif args.fix:
            print(indent, m.group(2))  # inserts single space between args
        elif indent + " " != m.group(1):
            print(f"{fi.filename()}:{fi.filelineno()} has indent {m.group(1)!r} but should have {indent + ' '!r}",
                  file=sys.stderr)
            broken = True

        if "*/" in line:
            indent = None

    else:
        if args.fix:
            print(line, end="")
        m = re.match("^([ \t]*)/\\*\\*", line)
        if not m or "*/" in line: continue
        indent = m.group(1)

if broken:
    sys.exit(1)

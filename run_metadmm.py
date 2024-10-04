
#!/usr/bin/env /Users/swati/opt/anaconda3/bin/python

import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))


from metaDMM.cli import cli

if __name__ == '__main__':
    cli()

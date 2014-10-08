#!/usr/bin/env python
#dump all metadata with given state to individual uuid names directories

import os
import sys
import CGHWSI

def main():
        if len(sys.argv) < 2:
		sys.stderr.write("must submit a CGHub analysis_id (uuid)\n")
		sys.exit(-1)
	uuid=sys.argv[1]
	sys.stdout.write("Processing uuid: %s\n" % uuid)
	data=CGHWSI.retrieve_analysis_attributes_for_uuid(uuid)
	error=CGHWSI.split_analysis_attributes(data,uuid)

if __name__ == '__main__':
	main()

import os

"""\
 We use this app to search for submited jobs at the farm which went on hold
 and consequently didn't output their root files. 
 This happens without errors, so it doesn't mean that the jobs failed
 simply the alocated resources were not sufficient to run them
 For now specific for electron scattering
 Authors: 
        Alon Sportes (Tel Aviv University) 
        Julia Tena Vidal (Tel Aviv University) 
Copyright:
   Copyright (c) 2003-2025, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org
"""

import os, optparse, glob, tarfile
import sys

def find_unmatched_scripts(jobs_dir,output_txt):
    PrinOut = False

    # Create file where we will list the missing files
    list_file = jobs_dir+"/list.txt"
    if os.path.exists(list_file) :
        os.remove(list_file)
    os.system("ls -1 "+jobs_dir+"/master-routine_validation_01-eScattering/ >" + list_file)

    # Obtain prefix
    with open(jobs_dir+"/grid_submission.xml", 'r') as xml_file:
        for xml_line in xml_file:
            if "<parallel>" not in xml_line:
                prefix_index = xml_line.rfind('/')
                prefix = xml_line[:prefix_index] + "/"
                break

    xml_file.close()

    # Initialize the lists
    gst_root_files = []
    sh_files = []

    # Read the input .txt file
    with open(list_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.endswith('.gst.root') or line.endswith('.xml'):
            # if line.endswith('.gst.root'):
                gst_root_files.append(line)
            elif line.endswith('.sh'):
                sh_files.append(line)

    print("Number of '.gst.root' files: " + str(len(gst_root_files)))
    print("Number of '.sh' files: " + str(len(sh_files)))

    print("\n")

    # List to store unmatched .sh files
    unmatched_sh_files = []

    # Check for unmatched .sh files
    for sh_file in sh_files:
        matched = False
        sh_file_base_name = os.path.splitext(sh_file)[0]

        if PrinOut:
            print(sh_file_base_name)

        for gst_file in gst_root_files:
            gst_file_base_name0 = os.path.splitext(gst_file)[0]
            gst_file_base_name = os.path.splitext(gst_file_base_name0)[0]

            if PrinOut:
                print(gst_file_base_name)

            if sh_file_base_name == gst_file_base_name:
                matched = True
                break
        if not matched:
            unmatched_sh_files.append(sh_file)

    # Write unmatched .sh files to the output .txt file
    if os.path.exists(output_txt) :
        os.remove(output_txt)
    with open(output_txt, 'w') as file:
        file.write('<parallel>\n')
        for sh_file in unmatched_sh_files:
            file.write(prefix + sh_file + '\n')
        file.write('</parallel>\n')

    print("Expected Number of '.gst.root' and '.sh' files: " + str(len(sh_files) * 2))
    print("Number of unmatched '.sh' files: " + str(len(unmatched_sh_files)))


if __name__ == "__main__":
    op = optparse.OptionParser(usage=__doc__)
    op.add_option("--jobs-topdir", dest="JOBSTD", default=os.getenv('PWD'), help="Top level dir for the job files (default: %default)")
    opts, args = op.parse_args()

    output_txt = opts.JOBSTD+"/grid_submission_unfinished.xml"  # output .xml file
    find_unmatched_scripts(opts.JOBSTD,output_txt)

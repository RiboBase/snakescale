import argparse
import os
import yaml
import subprocess
import glob
import shutil

def _download_references(ref_dict, target_folder):
    """
    Helper function for the actual download function
    """

    os.makedirs(target_folder, exist_ok = True)
    os.makedirs(target_folder + "/transcriptome", exist_ok = True)
    os.makedirs(target_folder + "/filter", exist_ok = True)

    for organism, contents in ref_dict.items():
        print("Downloading {} from {}".format(organism, contents["download_url"]) )
        download_file = contents["folder_name"] + ".tar.gz"

        if os.path.isfile(download_file):
            print("{} exists. Skipping the download.".format(download_file))
        else:
            download_res = subprocess.run(['curl', '-L', "--output", download_file, contents["download_url"]], 
                                          stdout = subprocess.PIPE)
            if download_res.returncode:
                raise Error("Error in downloading " + contents["download_url"])

        untar_res = subprocess.run(['tar', '-xzvf', download_file], 
                                      stdout = subprocess.PIPE)

        if untar_res.returncode:
            raise Error("Error in opening the archive " + download_file)

        # replace space with -
        # we use this to find the untar folder name
        untar_folder_name_start = "reference_" + organism.replace(" ", "-")

        folder_search_list = glob.glob("./" + untar_folder_name_start + "*")
        
        # There needs to be only one folder
        assert len(folder_search_list) == 1

        reference_untarred_folder = folder_search_list[0]

        for f in ["filter", "transcriptome"]:
            f_contents = glob.glob(reference_untarred_folder + "/" + f + "/*")
            for x in f_contents:
                shutil.move(x, target_folder + "/" + f + "/")

        shutil.rmtree(reference_untarred_folder)


def download_reference(reference_yaml, target_folder):
    print("Downloading the references....")

    #which_proc = subprocess.Popen(['which', 'curlerr'])
    #outs, errs = which_proc.communicate()
    #print("which result:" ,errs)
    # args, returncode, stdout, stderr
    whic_res = subprocess.run(['which', 'curl'], stdout = subprocess.PIPE)

    if whic_res.returncode: 
        raise FileNotFoundError("Could not find the executable curl")

    if not os.path.isfile(reference_yaml):
        raise FileNotFoundError("Could not find the yaml file:", reference_yaml)

    with open( reference_yaml, 'rt' ) as yaml_stream:
        yaml_contents = yaml.load(yaml_stream, Loader=yaml.FullLoader)

    organisms = list(yaml_contents.keys())
    print("The following organisms are going to be downloaded:")
    print("  -", "\n  - ".join(organisms))

    _download_references(yaml_contents, target_folder)

    #arrange_references(yaml_contents)


def get_parameters():
    parser = argparse.ArgumentParser(
                  formatter_class = argparse.RawTextHelpFormatter,
                  description = \
                      
    """
    Download & Organize RiboFlow Refrences
    The references are read from the yaml file. They are gathered under the folder given in --target
    """)

    parser.add_argument("--target",
                        type     = str,
                        required = False,
                        default  = "./reference",
                        help     = "path of the folder to download references to")

    parser.add_argument("--yaml",
                        type     = str,
                        required = False,
                        default  = "references.yaml",
                        help     = "path of the folder to download references to")
    

    args = parser.parse_args()
    return args


def main():
    params  = get_parameters()
    download_reference(target_folder = params.target,
                       reference_yaml = params.yaml)
    print("References have been gathered in the folder:", params.target)

if __name__ == "__main__":
    main()


# Sample Yaml File
# ###########################################################
# "homo sapiens":
#   download_url: "https://github.com/RiboBase/reference_homo-sapiens/archive/refs/tags/v1.0.tar.gz"
#   official_name: "homo sapiens"
#   folder_name: "human"
#   filter: "filter/human/human_rtRNA*"
#   transcriptome: "transcriptome/human/appris_human_v2_selected*"
#   regions: "transcriptome/human/appris_human_v2_actual_regions.bed"
#   transcript_lengths: "transcriptome/human/appris_human_v2_transcript_lengths.tsv"

# "mus musculus":
#   download_url: "https://github.com/RiboBase/reference_mus-musculus/archive/refs/tags/v1.0.tar.gz"
#   official_name: "mus musculus"
#   folder_name: "mouse"
#   filter: 'filter/mouse/mouse_rtRNA*'
#   regions: 'transcriptome/mouse/appris_mouse_v2_actual_regions.bed'
#   transcript_lengths: 'transcriptome/mouse/appris_mouse_v2_transcript_lengths.tsv'
#   transcriptome: 'transcriptome/mouse/appris_mouse_v2_selected*'    
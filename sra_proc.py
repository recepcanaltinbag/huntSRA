
from logging import exception
from urllib import response
from Bio import Entrez
import pandas as pd
import subprocess
import shutil
import os
import re
import time
import glob
import random
import requests

from tqdm import tqdm
import urllib.request

Entrez.email = "lalala@gmail.com"


class DownloadProgressBar(tqdm):
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


def download_url(url, output_path):
    with DownloadProgressBar(unit='B', unit_scale=True,
                             miniters=1, desc=url.split('/')[-1]) as t:
        urllib.request.urlretrieve(url, filename=output_path, reporthook=t.update_to)


def SRA_id_to_ACC(sra_id):
    print('Find ACC: ', sra_id)
    esummary_handle = Entrez.esummary(db="sra", id=sra_id, report="full")
    esummary_record = Entrez.read(esummary_handle, validate=False)
    esummary_handle.close()
    print('end search')
    start_c = 'Run acc="'
    end_c = '" total_spots='
    SRR_ACC = re.search(re.escape(start_c)+"(.*)"+re.escape(end_c),esummary_record[0]['Runs']).group(1)
    return SRR_ACC


#threshold size is Mb
def is_content_size_lower(url, threshold_size, lower_threshold_size):
    response = requests.head(url, allow_redirects=True)
    if 'Content-Length' in response.headers:
        #return how many Mb
        the_size = int(response.headers['Content-Length'])/(1024*1024)
        if the_size > threshold_size or the_size < lower_threshold_size:
            print('file size is ', the_size ,' is bigger than threshold: ', threshold_size)
            return False
        else:
            return True
    else:
        raise Exception("No Content Length in request")


def run_sra_down2(sra_acc, out_path, threshold_size, lower_threshold_size):

    if not os.path.exists("sra_download"):
        os.mkdir("sra_download")
    if not os.path.exists(f"sra_download/{sra_acc}"):
        os.mkdir(f"sra_download/{sra_acc}")

    the_url = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{sra_acc}/{sra_acc}"
    #args2 = f"fastq-dump --split-3 {out_path}/*.sra --outdir {out_path}"
    args2 = f"parallel-fastq-dump --sra-id {sra_acc} --split-3 --threads 24"
    print(out_path)

    if is_content_size_lower(the_url, threshold_size, lower_threshold_size) == False:
        return False

    print('Content Size is reasonable')
    
    print(out_path + "/" + sra_acc + ".sra")
    if not os.path.exists(out_path + "/" + sra_acc + ".sra"):
        print('downloading')
        download_url(the_url, out_path + "/" + sra_acc + ".sra")
    if not os.path.exists(out_path + "/" + sra_acc + "_1.fastq"):   
        print('fastqing')
        cwd = os.getcwd()
        os.chdir(cwd+f'/{out_path}')
        try:
            my_process = subprocess.run(args2, shell=True, text=True, capture_output=True)
            os.chdir(cwd)
        except Exception as e:
            my_process = 0
            os.chdir(cwd)
            print('url download error: ', e)
            return False

    print('Start for Fastq Dump')
    #my_process = subprocess.run(args2, shell=True, text=True, capture_output=True)

    print('End of Fastq Dump')
    #if os.path.exists(f"{out_path}/*.sra"):
    #    shutil.rmtree(f"{out_path}/*.sra")
    #    print('SRA Files deleted.')
    return True


def run_sra_down(sra_acc, out_path):
    args = f"sra-downloader {sra_acc} --save-dir {out_path}"
    my_process = subprocess.run(args, shell=True, text=True, capture_output=True)
    if my_process.returncode != 0:
        raise Exception("Error in SRA Download")
    return True


def control_fastq_files(file_list):
    if len(file_list) > 2:
        new_list = []
        for file_name in file_list:
            if '_1' in file_name:
                new_list.append(file_name)
            if '_2' in file_name:
                new_list.append(file_name)
        return new_list
    else:
        return file_list


def run_sam_pipeline(sra_file_path, query_file, q_name):

    if not os.path.exists("sra_cov"):
        os.mkdir("sra_cov")
    
    if not os.path.exists("sra_temp"):
        os.mkdir("sra_temp")
    

    fastq_files = glob.glob(sra_file_path + "/*.fastq")
    #stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    print("fastq files: ", fastq_files)

    if len(fastq_files) == 0:
        raise Exception("No fastq files are found, skipping!!")

    fastq_files = control_fastq_files(fastq_files)
    fastq_files = (" ").join(fastq_files)

    args = f"bwa index -b 100000000 {query_file}"
    my_process = subprocess.run(args, shell=True, text=True, capture_output=True)
    if my_process.returncode != 0:
        raise Exception("Error in bwa Indexing", my_process)
    else:
        print('BWA Indexing is successful')
    args = f"bwa mem -t 24 {query_file} {fastq_files} > {sra_file_path}/aln{q_name}.sam" 
    my_process = subprocess.run(args, shell=True, text=True, capture_output=True)
    if my_process.returncode != 0:
        raise Exception("Error in bwa Mem", my_process)
    else:
        print('BWA Mem is successful')
    #    if os.path.exists('sra_download'):
    #        shutil.rmtree('sra_download')
    #        print('SRA Download Folder Deleted!')

    args = f"samtools view --threads 24 -S -b {sra_file_path}/aln{q_name}.sam > {sra_file_path}/bamfile{q_name}.bam" 
    my_process = subprocess.run(args, shell=True, text=True, capture_output=True)
    if my_process.returncode != 0:
        raise Exception("Error in samtools view", my_process)
    else:
        print('Samtools View is successful')
    args = f"samtools sort --threads 24 -m 800M {sra_file_path}/bamfile{q_name}.bam -o {sra_file_path}/sorted{q_name}.bam"
    my_process = subprocess.run(args, shell=True, text=True, capture_output=True)
    if my_process.returncode != 0:
        raise Exception("Error in samtools sort", my_process)
    else:
        print('Samtools sort is successful')
    args = f"samtools coverage {sra_file_path}/sorted{q_name}.bam -o {sra_file_path}/cov_info{q_name}.txt"
    my_process = subprocess.run(args, shell=True, text=True, capture_output=True)
    if my_process.returncode != 0:
        raise Exception("Error in Samtools Coverage", my_process)
    else:
        print('Samtools Coverage is successful')
    
    return True

def read_cov_file(file_name, threshold_cov=90.0):
    try:
        cov_result = pd.read_table(file_name)
    except pd.errors.EmptyDataError:
        print(" is empty and has been skipped. %s" % cov_result)
        return False
    else:
        print(cov_result)
        df_filtered = cov_result[(cov_result['coverage'] >= threshold_cov)]
        if len(df_filtered) > 0:
            print('mapping result true')
            return True
        else:
            print('mapping result false')
            return False
    return False

def clean_temp_files(out_path):
    try:
        if os.path.exists("sra_cov"):
            shutil.rmtree("sra_cov")
        if os.path.exists("sra_temp"):
            shutil.rmtree('sra_temp')
        if os.path.exists(out_path):
            shutil.rmtree(out_path)
        if os.path.exists('sra_download'):
            shutil.rmtree('sra_download')
    except:
        print('DELETE FILE ERROR')

def control_programs():
    my_process = subprocess.run("samtools help", shell=True, text=True, capture_output=True)
    if my_process.returncode != 0:
        raise Exception("Samtools could not found, build and copy to PATH")
    my_process = subprocess.run("bwa index sra_query/QxyA_IS1380.fasta", shell=True, text=True, capture_output=True)
    my_process = subprocess.run("bwa index sra_query/QxyA_Tn3.fasta", shell=True, text=True, capture_output=True)

    if my_process.returncode != 0:
        raise Exception("BWA could not found sssssss, copy bam to PATH")
    my_process = subprocess.run("fasterq-dump -h", shell=True, text=True, capture_output=True)
    if my_process.returncode != 0:
        raise Exception("SRA Toolkit couldnt found, install sra-toolkit to the PATH")
    my_process = subprocess.run("sra-downloader -h", shell=True, text=True, capture_output=True)
    if my_process.returncode != 0:
        raise Exception("sra-downloader couldnt found, pip3 install sra-downloader")
    return True


    

def main():

    #Control Programs
    try:
        control_programs()
    except Exception as e:
        print('Some programs are missing:\n', e)
    else:
        print('all programs are installed, continuing')
        threshold_size = 50000 #if file is bigger than 50000 Mb, dont donwload
        threshold_size = 10000 #if file is bigger than 50000 Mb, dont donwload

        lower_threshold_size = 25 #if tyhe file is smaller than 10 Mb, dont, download
        file_path = "SRA_ids_qxyA_Tn3_IS1380.tab"  # Dosya yolunu buraya yazın

        sra_ids = []
        with open(file_path, "r", encoding="utf-8") as f:
            sra_ids = [line.strip() for line in f]
        print(len(sra_ids))

        query_Tn3_file = "sra_query/QxyA_Tn3.fasta"
        query_IS1380_file = "sra_query/QxyA_IS1380.fasta"
        out_dir = "sra_download"
        sra_counter = 0
        for sra_id in sra_ids:
            sra_counter = sra_counter + 1
            print("Start Analysis for: ", sra_id, "len of SRAs:", sra_counter,'/',len(sra_ids))
            named_tuple = time.localtime() # get struct_time
            time_string = time.strftime("%m/%d/%Y, %H:%M:%S", named_tuple)
            print(time_string)

            sra_acc = sra_id    
            #sra_acc = SRA_id_to_ACC(sra_id)
            print(sra_acc)
            out_path = os.path.join(out_dir, sra_acc)
            #clean_temp_files(out_path)
            try:
                run_sra_down2(sra_acc, out_path, threshold_size, lower_threshold_size)
            except Exception as e:
                print("SRA DOWNLOAD ERROR: ", sra_id)
                print('ERROR:', e)
                #clean_temp_files(out_path)
            else:
                print('Download completed, continue to analysis')
                if not os.path.exists(f"{out_path}/cov_infoTn3.txt"):
                    try:
                        run_sam_pipeline(out_path, query_Tn3_file, "Tn3")
                        run_sam_pipeline(out_path, query_IS1380_file, "IS1380")
                    except Exception as e:
                        print('ERROR:', e)
                        print("SRA ANALYSIS ERROR: ", sra_id)
                    else:
                        if read_cov_file(f"{out_path}/cov_infoTn3.txt"):
                            bio_status = "LOOKED_FOUND Tn3"
                            input()
                            #Bio_DB.insert_SRA_Info(sra_id, project_id, bio_status)
                            #Bio_DB.insert_SRA_Found(sra_id, project_id)              
                            print("*****-> Found and Added DB")
                        if read_cov_file(f"{out_path}/cov_infoIS1380.txt"):
                            bio_status = "LOOKED_FOUND IS1380"
                            input()
                            #Bio_DB.insert_SRA_Info(sra_id, project_id, bio_status)
                            #Bio_DB.insert_SRA_Found(sra_id, project_id)              
                            print("*****-> Found and Added DB")

            #clean_temp_files(out_path)
            print('End of analysis')


if __name__ == '__main__':
    while True:
        try:
            main()
        except:
            print('An error occured, will try again')
            time.sleep(30)
else:
    print('The main() function did not execute')

#read sam coverage output file and parse according to the result, threshold is 90 as default
#print(read_cov_file("sra-test/cov_info.txt"))
#print(re.search(re.escape(start)+"(.*)"+re.escape(end),s).group(1))
#sra-downloader SRR4235047 --save-dir sra-test
#my_process = subprocess.run(args, shell=True, executable='/bin/bash', text=True, check=True,  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
# ID to RUNs
#sra-downloader SRR11362851 --save-dir sra-test
#bwa mem OxyBAC.fasta *1.fastq.gz *2.fastq.gz > aln.sam
# It needed sorted bam file to calculate coverage 
#samtools coverage sorted.bam -o cov_info.txt
#handle = Entrez.efetch(db="nucleotide", id="AY851612", rettype="gb", retmode="text")
#esummary_handle = Entrez.esummary(db="bioproject", id=my_id, report="full", retmode="xml")
#esummary_record = Entrez.read(esummary_handle, validate=False)
#esummary_handle.close()
#print(esummary_record['DocumentSummarySet']['DocumentSummary'][0])
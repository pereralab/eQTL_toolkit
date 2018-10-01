#!/usr/bin/python3

import subprocess
import argparse
import get_eqtl
import os

def cd(cd_path):
        saved_path = os.getcwd()
        os.chdir(cd_path)
        yield
        os.chdir(saved_path)


def runCis(vcf, bed, cov, dirs, prefix, analysis):
        for i in range(1,30):
                print("cis-eQTL of chunck " + str(i)+ "/30")
                cmd = "cp " + dirs + "/src/submit.sh "+ dirs +"/src/" + prefix + "_cis_eQTL_"+ analysis + "_" + str(i)+".sh"
                
                os.system(cmd)
                cmd1 = "cd " + dirs + "/eQTL " + "\n\n\n\n"

                
                #cmd2 = "QTLtools_1.1_Ubuntu12.04_x86_64 cis \\\n" \
                cmd2 = "/projects/b1047/zhong/software/FastQTL/bin/fastQTL.static  \\\n" \
                + "--vcf " + dirs + "/data/" + vcf + " \\\n" \
                + "--bed " + dirs + "/data/" + bed + " \\\n" \
                + "--cov " + dirs + "/data/" + cov + " \\\n" \
                + "--chunk " + str(i) + " 30" + " \\\n" \
                
                
                if analysis == "nominal":
                        cmd2 = cmd2  \
                        + "--out " + prefix + "_cis_" + str(i) + "_30.txt" \
                        #+ "--nominal 1 \\\n"

                elif analysis == "permute":
                        cmd2 = cmd2 + "--permute 1000 \\\n" \
                        + "--out " + prefix + "_permute_" + str(i) + "_30.txt"

                with open(dirs +"/src/" + prefix + "_cis_eQTL_"+ analysis + "_" + str(i)+".sh" , "a") as myfile:
                        myfile.write(cmd1)
                        myfile.write(cmd2)
                myfile.close()
                os.system("chmod +x " + dirs +"/src/" + prefix + "_cis_eQTL_"+ analysis + "_"+str(i)+".sh")
                os.system("msub " + dirs + "/src/" + prefix + "_cis_eQTL_"+ analysis + "_" + str(i)+".sh")

if __name__=='__main__':
        parser = argparse.ArgumentParser(description='eQTL analysis with QTLTools')
        parser.add_argument('vcf', help='vcf file')
        parser.add_argument('bed', help='bed file')
        parser.add_argument('cov', help='covariate file')
        parser.add_argument('--working_dir', help='directory', default='./')
        parser.add_argument('--analysis',default=None, help='analysis')
        parser.add_argument('--prefix', default='eQTL', help='result file prefix')
        args = parser.parse_args()


        vcf = args.vcf
        bed = args.bed
        cov = args.cov
        prefix = args.prefix
        dirs = args.working_dir
        
        if args.analysis == "nominal":
                
                print("cis-eQTL mapping nomianl mode")
                runCis(vcf, bed, cov, dirs, prefix, "nominal")


        elif args.analysis == "permute":
                print("cis-eQTL mapping permutation mode")
                runCis(vcf, bed, cov, dirs, prefix, "permute")

        elif args.analysis == "call":
                print("call significant eQTLs")
                
                cd(dirs + "/eQTL/")

                cmd = "cat " + prefix + "_cis_*_30.txt | gzip -c > " + prefix + "_cis_full.txt.gz"
                #os.system(cmd
                
                cmd = "cat " + prefix + "_permute_*_30.txt | gzip -c > " +  prefix + "_permute_full.txt.gz"
               # os.system(cmd)
                
                #cmd = "Rscript " + dirs + "/src/runFDR_cis.R " + prefix + "_permute_full.txt.gz" + " 0.05 " + prefix
                cmd = "Rscript " + dirs + "/src/calulateNominalPvalueThresholds.R " + prefix + "_permute_full.txt.gz" + " 0.05 " + prefix
                #os.system(cmd)

                egene = get_eqtl.threshold(prefix)
                #print(egene)
                cmd = "awk \'NR==FNR {a[$1];next} ($1 in a) {print $0}\' " \
                + prefix + ".significant.txt " \
                + "<(zcat " + prefix + "_cis_full.txt.gz)" \
                + "> " + prefix + "_cis_full.sig.txt"
                print(cmd)
                #subprocess.call(cmd, shell=True)

                get_eqtl.call(egene, prefix + "_cis_full.sig.txt", prefix + ".significant_eQTLs.txt") 


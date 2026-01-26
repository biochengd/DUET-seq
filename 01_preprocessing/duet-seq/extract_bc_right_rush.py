#!/s1/Liuyang/biosoft/miniconda3/bin/python
# coding: utf-8
import pyfastx
from fuzzysearch import find_near_matches
import os
import argparse

def read_txt_file(txt_file_path):
    with open(txt_file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]

# 序列比对
def sequence_alignment(txt_lines, fastq_sequence, max_l_dist=0):
    for barcode in txt_lines:
        match = find_near_matches(barcode, fastq_sequence[13:24], max_l_dist=0)
        if match:
            barcode3_end = match[0].end+13
            return barcode3_end
            break

# 运行主程序
# 主要流程
def extract_bc_and_save(r1_fastq_path, r2_fastq_path, folder_path):
    txt_file_path="/s1/SHARE/barcode3.txt"
    txt_lines = read_txt_file(txt_file_path)
    r1_sequences =  [record for record in pyfastx.Fastq(r1_fastq_path)]
    r2_sequences = [record for record in pyfastx.Fastq(r2_fastq_path)]

    # 判断文件夹是否存在
    if not os.path.exists(folder_path):
        # 如果文件夹不存在，创建新的文件夹
        os.makedirs(folder_path)

    # 定义output的名称
    Barcode = './'+folder_path+'/barcode.fq'
    R1_filter = './'+folder_path+'/R1_filter.fq'

    barcode_fq = open(f'{Barcode}', 'w')
    r1_filter_fq = open(f'{R1_filter}', 'w')

    for r1_seq, r2_seq in zip(r1_sequences, r2_sequences):
        barcode3_end = sequence_alignment(txt_lines, r2_seq.seq)
        barcode_1 = ""
        barcode = ""
        if barcode3_end:
                # 简化的条形码处理逻辑
                # 仅示例，具体逻辑可能需要根据实际情况调整
            barcode_1 = r2_seq.seq[:3]  # 默认情况
            special_cases = {
                22: {"GATA": "GAA"},
                23: {"ATTGA": "TTG"}
            }
                # 检查特殊情况
            for key, value in special_cases.get(barcode3_end, {}).items():
                if r2_seq.seq.startswith(key):
                    barcode_1 = value
                    break
            if barcode3_end == 24:
                barcode_1 = r2_seq.seq[3:6]
                ### barcode3 取后六位 barcode2 取前7位
            barcode = f"{barcode_1}{r2_seq.seq[barcode3_end-17:barcode3_end-10]}{r2_seq.seq[barcode3_end-6:barcode3_end]}"
            umi = r2_seq.seq[barcode3_end:barcode3_end+10]

            # 写入对应文件
            barcode_fq.write(f"@{r2_seq.name}\n{barcode}GG{umi}\n+\n{r2_seq.qual[:5]}{r2_seq.qual[barcode3_end-17:barcode3_end-10]}{r2_seq.qual[barcode3_end-6:barcode3_end+10]}\n")
            r1_filter_fq.write(f"@{r1_seq.name}\n{r1_seq.seq[:90]}\n+\n{r1_seq.qual[:90]}\n")

    # 保存新的 FASTQ 文件
    # 关闭文件
    barcode_fq.close()
    r1_filter_fq.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract barcodes and save to files.")
    parser.add_argument('-r1', '--r1_fastq_path', required=True, help='Path to R1 fastq file')
    parser.add_argument('-r2', '--r2_fastq_path', required=True, help='Path to R2 fastq file')
    parser.add_argument('-o', '--output_folder', required=True, help='Output folder path')

    args = parser.parse_args()
    extract_bc_and_save(args.r1_fastq_path, args.r2_fastq_path, args.output_folder)

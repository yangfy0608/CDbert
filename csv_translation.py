import csv
import os
import sys
import argparse
from datetime import datetime
import random
csv.field_size_limit(100000000)

parser = argparse.ArgumentParser(
    description='Filter CSV file by specified date')

parser.add_argument('-i', '--input', type=str,
                    help='Input Folder name', 
                    required=True)

parser.add_argument('-t', '--time', type=str,
                    help='Input time with "DD-MON-YYYY", like "15-MAY-2016"', 
                    required=True)

parser.add_argument('-min', '--minlength', type=int,
                    help='the min length of the gene', 
                    required=False, 
                    default=150)

parser.add_argument('-max', '--maxlength', type=int,
                    help='the length of the gene u want', 
                    required=False, 
                    default=251)

parser.add_argument('-n', '--numseqs', type=int,
                    help='number of sequences to randomly select', 
                    required=False, 
                    default=None)

args = parser.parse_args()

# 指定日期
specified_date_object = datetime.strptime(args.time, "%d-%b-%Y")
specified_timestamp = specified_date_object.timestamp()

# 读取CSV文件
for filename in os.listdir(args.input):
    if filename.endswith('.csv'):
        filepath = os.path.join(args.input, filename)
        with open(filepath, "r") as csvfile,                                       \
        open(filepath.replace(".csv", "_vibe_1.fastq"), "w", newline='') as fastqfile1, \
        open(filepath.replace(".csv", "_vibe_2.fastq"), "w", newline='') as fastqfile2, \
        open(filepath.replace(".csv", "_dvf.fasta"), "w", newline='') as fastafile:
            num_lines = 0
            reader = csv.DictReader(csvfile)
            seqs = []
            for row in reader:
                # 将日期字符串转换为日期对象
                date_string = row['submit date']
                # print(row)
                date_object = datetime.strptime(date_string, "%d-%b-%Y")
                # 将日期对象转换为Unix时间戳
                timestamp = date_object.timestamp()
                # 比较Unix时间戳并输出符合条件的行
                if timestamp > specified_timestamp and len(row['sequence']) >= 2*args.minlength:
                    # print(date_string)
                    seqs.append(row)
                else:
                    pass

            # 随机选择指定数量的行（如果有）
            if args.numseqs is  None:
                print("需要给定抽取比例(0,100)")
                sys.exit(1)
            elif int(args.numseqs * len(seqs) / 100) < 1:
                print("在",filename,"中满足条件的序列少于需求")
                sys.exit(1)
            else:
                pre_num = int(args.numseqs * len(seqs) / 100)
                seqs = random.sample(seqs, pre_num)

            for i, row in enumerate(seqs):
                seq_id = row['name']
                seq_data = row['sequence']
                seq_length = len(seq_data)

                if seq_length < args.maxlength*2:
                    split_idx = seq_length // 2
                    seq2 = seq_data[:split_idx] if split_idx >= args.maxlength else seq_data[:args.maxlength]
                    seq3 = seq_data[split_idx:split_idx+args.maxlength] if seq_length - \
                        split_idx >= args.maxlength else seq_data[split_idx:]
                    qual_score2 = '~' * len(seq2)
                    qual_score3 = '~' * len(seq3)

                    # Write FASTQ format
                    fastqfile1.write(
                        '@{}\n{}\n+\n{}\n'.format(seq_id, seq2, qual_score2))
                    fastqfile2.write(
                        '@{}\n{}\n+\n{}\n'.format(seq_id, seq3, qual_score3))
                    fastafile.write(">" + seq_id + "\n" + seq2+seq3 + "\n")

                else:
                    # 截取sequence
                    split_idx = random.randint(0, seq_length - args.maxlength*2)
                    seq2 = seq_data[split_idx:split_idx+args.maxlength] if seq_length - \
                        split_idx >= args.maxlength else seq_data[split_idx:]
                    seq3 = seq_data[split_idx+args.maxlength:split_idx+args.maxlength*2] if seq_length - \
                        split_idx-args.maxlength >= args.maxlength else seq_data[split_idx+args.maxlength:]

                    # 生成相应长度的质量分数
                    qual_score2 = '~' * len(seq2)
                    qual_score3 = '~' * len(seq3)

                    # Write FASTQ format
                    fastqfile1.write(
                        '@{}\n{}\n+\n{}\n'.format(seq_id, seq2, qual_score2))
                    fastqfile2.write(
                        '@{}\n{}\n+\n{}\n'.format(seq_id, seq3, qual_score3))
                    fastafile.write(">" + seq_id + "\n" + seq2+seq3 + "\n")
            print("成功从",filename,"提取", pre_num, "条碱基对序列并转化为适用于vibe_predict的fastq文件")
            print("成功从",filename,"提取", pre_num, "条碱基对序列并转化为适用于dvf_predict的fasta文件")

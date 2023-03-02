import pysam
import argparse
from pyfaidx import Fasta
import time
import copy

parser = argparse.ArgumentParser()
parser.add_argument("hg19_bam", help='bam file')
parser.add_argument("fasta", help='fasta file')
args = parser.parse_args()
localtime = time.asctime(time.localtime(time.time()))
print ('Start',localtime)
merge_distance=1000
low_rate=0.8
high_rate=1.2
mini_SV_LEN=40

main_chr=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
          'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrY','chrX',
          '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12','13', '14', '15', '16', '17',
          '18', '19', '20', '21', '22', 'Y', 'X']

samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
DEL_final_res=[]
#point= open('site2.txt', 'w')
record=[]
#detect the candidate DEL
for r in samfile.fetch():
  if r.mapq>20 and r.is_supplementary==False:
      chr = r.reference_name
      start_pos=r.reference_start
      cigar=r.cigarstring
      read_str=r.query_sequence
      i=0
      pos=start_pos
      pos1=0
      pos2=0
      while i<len(cigar):
          num=''
          while cigar[i].isdigit():
              num=num+cigar[i]
              i+=1
          num=int(num)
          if cigar[i]=='M':
              pos=pos+num
          elif cigar[i]=='D':
              pos1=pos
              pos=pos+num
              pos2=pos
              if num>=30:
                 record.append([chr, pos1, pos2])
          i+=1
samfile.close()
#merge the candidate DEL
res=[]
intervals =list(sorted(record))
chr=intervals[0][0]
low = intervals[0][1]
high = intervals[0][2]
read_num=1
for i in range(1, len(intervals)):
    newSV_len=intervals[i][2]-intervals[i][1]
    if (high+merge_distance) >= intervals[i][1] and intervals[i][0]==chr:
        low = min(intervals[i][1],low)
        high = max(intervals[i][2],high)
        read_num += 1
    else:
        res.append([chr,low, high,read_num])
        chr=intervals[i][0]
        low = intervals[i][1]
        high = intervals[i][2]
        read_num=1
res.append([chr,low, high,read_num])
newres=[]
for site in res:
    if int(site[3])>=3 and (site[2]-site[1])>=30:
        newres.append([site[0],site[1],site[2]])
res_region=newres

newDEL=[]
samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
for region in res_region:
    chr = region[0]
    sv_pos1 = region[1]
    sv_pos2 = region[2]
    read_del=[]
    zero_DEL_num = 0
    side_len = 1000
    coverage=0
    for r in samfile.fetch(chr, sv_pos1, sv_pos2):
        if r.mapq > 20 and r.is_supplementary == False and \
                r.reference_start<(sv_pos1 - 1100) and r.reference_end>(sv_pos2 + 1100):
            coverage+=1
            range1=sv_pos1 - side_len
            range2 = sv_pos2 + side_len
            DEL_num=0
            DEL_len=0
            DEL_pos=[]
            start_pos = r.reference_start
            cigar = r.cigarstring
            read_str = r.query_sequence
            i = 0
            pos = start_pos
            pos1 = 0
            pos2 = 0
            while i < len(cigar):
                num = ''
                while cigar[i].isdigit():
                    num = num + cigar[i]
                    i += 1
                num = int(num)
                if cigar[i] == 'M':
                    pos = pos + num
                elif cigar[i] == 'D':
                    pos1 = pos
                    pos = pos + num
                    pos2 = pos
                    if pos1>range1 and pos2<range2:
                        if num>=40:
                            DEL_pos.append([chr, pos1, pos2])
                            DEL_num += 1
                            DEL_len+=num
                i += 1
            if DEL_num!=0:
                read_del.append([DEL_num, DEL_len, DEL_pos])
            else:
                zero_DEL_num+=1
            if coverage>150:
                break

    if len(read_del)==0 or float(len(read_del))/(float(len(read_del))+float(zero_DEL_num))<0.2:
        continue
    read_del = sorted(read_del, key=lambda x: x[1])
    noise_len = int(len(read_del) / 20)
    if noise_len != 0:
        read_del = read_del[noise_len:-noise_len]

    sv_len_num = []
    sv_len = read_del[0][1]
    read=read_del[0]
    num = 1
    i = 1
    while i < (len(read_del)):
        if read_del[i][1] == sv_len:
            if read_del[i][0]<read[0]:
                read=read_del[i]
            num += 1
        else:
            sv_len_num.append([sv_len,num,read])
            sv_len = read_del[i][1]
            read=read_del[i]
            num = 1
        i += 1
    sv_len_num.append([sv_len, num, read])

    cengci_cluster = []
    if len(sv_len_num) <= 2:
        cengci_cluster.append(sv_len_num)
    while len(sv_len_num) > 2:
        distance = []
        distance_rate=[]
        i = 1
        while i < len(sv_len_num):
            node_dis = sv_len_num[i][0] - sv_len_num[i - 1][0]
            node_dis_rate = float(sv_len_num[i][0] - sv_len_num[i - 1][0])/float(sv_len_num[i - 1][0])
            distance.append(node_dis)
            distance_rate.append(node_dis_rate)
            i += 1

        i = 1
        min_dis = distance[0]
        min_dis_rate=distance_rate[0]
        sign = 0
        while i < len(distance):
            if distance[i] < min_dis:
                min_dis = distance[i]
                min_dis_rate=distance_rate[i]
                sign = i
            i += 1

        if min_dis_rate<=0.2 or min_dis<=200:
            node_support_num = sv_len_num[sign][1] + sv_len_num[sign + 1][1]
            node_len = (float(sv_len_num[sign][0]) * sv_len_num[sign][1] + float(sv_len_num[sign + 1][0]) *
                        sv_len_num[sign + 1][1])
            node_len = node_len / node_support_num
            if sv_len_num[sign][2][0] < sv_len_num[sign + 1][2][0]:
                node_read = sv_len_num[sign][2]
            elif sv_len_num[sign][2][0] > sv_len_num[sign + 1][2][0]:
                node_read = sv_len_num[sign + 1][2]
            else:
                if abs(sv_len_num[sign][2][1] - node_len) < abs(sv_len_num[sign + 1][2][1] - node_len):
                    node_read = sv_len_num[sign][2]
                else:
                    node_read = sv_len_num[sign + 1][2]
            del sv_len_num[sign:sign + 2]
            sv_len_num.insert(sign, [node_len, node_support_num, node_read])
            sv_len_num_2 = copy.deepcopy(sv_len_num)
            cengci_cluster.append(sv_len_num_2)
        else:
            break

    if len(sv_len_num) == 1:
        chr = sv_len_num[0][2][2][0][0]
        del_pos1 = sv_len_num[0][2][2][0][1]
        del_pos2 = sv_len_num[0][2][2][-1][2]
        del_len = sv_len_num[0][0]
        newDEL.append([chr, del_pos1, del_pos2, sv_len_num[0][1], del_len])
    elif len(sv_len_num) == 2:
        cengci_cluster.pop()

        sv_len_num = sorted(sv_len_num, key=lambda x: x[1])
        sv_len_num.reverse()
        sv_len_num = sv_len_num[0:2]

        total_num = sv_len_num[0][1] + sv_len_num[1][1]
        support_read1=sv_len_num[0][1]
        support_read2=sv_len_num[1][1]
        rate1 = float(sv_len_num[0][1]) / float(total_num)
        rate2 = float(sv_len_num[1][1]) / float(total_num)
        del1=[sv_len_num[0][2][2][0][0],sv_len_num[0][2][2][0][1],sv_len_num[0][2][2][-1][2],sv_len_num[0][1],sv_len_num[0][0]]
        del2 = [sv_len_num[1][2][2][0][0], sv_len_num[1][2][2][0][1], sv_len_num[1][2][2][-1][2],sv_len_num[1][1],sv_len_num[1][0]]
        average_len = sv_len_num[0][1] * sv_len_num[0][0] + sv_len_num[1][1] * sv_len_num[1][0]
        average_len = float(average_len) / float(total_num)
        len_dif = abs(sv_len_num[1][0] - sv_len_num[0][0])
        while rate1<0.2 or rate2<0.2:
            if len(cengci_cluster)==0:
                rate1 = 0.5
                rate2 = 0.5
                len_dif = 0
                break
            cluster=sorted(cengci_cluster[-1], key=lambda x: x[1])
            total_num=cluster[-1][1]+cluster[-2][1]
            support_read1 = cluster[-1][1]
            support_read2 = cluster[-2][1]
            rate1=float(cluster[-1][1])/float(total_num)
            rate2 = float(cluster[-2][1]) / float(total_num)
            del1 = [cluster[-1][2][2][0][0], cluster[-1][2][2][0][1], cluster[-1][2][2][-1][2],cluster[-1][1],cluster[-1][0]]
            del2 = [cluster[-2][2][2][0][0], cluster[-2][2][2][0][1], cluster[-2][2][2][-1][2], cluster[-2][1],
                    cluster[-2][0]]
            average_len = cluster[-1][1] * cluster[-1][0] + cluster[-1][1] * cluster[-1][0]
            average_len = float(average_len) / float(total_num)
            len_dif = abs(cluster[-1][0] - cluster[-2][0])
            cengci_cluster.pop()
            if cengci_cluster==0:
                break

        len_dif_value = 20 + float(0.005 * average_len)
        if len_dif > len_dif_value:
            if rate1 >= 0.2 and  support_read1 >= 3:
                newDEL.append(del1)
            if rate2 >= 0.2 and  support_read2 >= 3:
                newDEL.append(del2)
        else:
            if (rate1+rate2)>=0.2 and (support_read1+support_read2)>=3:
                newDEL.append(del1)

samfile.close()
short_DEL=newDEL
SV_LEN=300
#detect long DEL
samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
long_DEL_record=[]
for r in samfile.fetch():
    if r.mapq > 20 and r.is_supplementary == False:
        sign = ""
        for tag in r.tags:
            if str(tag[0])=="SA":
                sign="SA"
                supple_align=tag
        # judge the two alignment
        if sign == "SA":
            # check primary alignment
            chr1 = r.reference_name
            cigar1 = r.cigarstring
            read = r.query_sequence
            other_align = supple_align[1].split(";")
            supple = other_align[0].split(",")
            i = 1
            while i < (len(other_align) - 1):
                m = other_align[i].split(",")
                if int(m[4]) > int(supple[4]):
                    supple = m
                i += 1
            chr2 = supple[0]
            # judge the primary alignment and supplementary alignment pos
            if r.reference_start <= int(supple[1]) and chr1 == chr2:
                pos1 = r.reference_end
                pos2 = int(supple[1])
                direct2 = supple[2]
                cigar2 = supple[3]
                s2 = cigar2.split("S")[0]
                if 'M' in s2:
                    s2 = 0
                else:
                    s2 = int(s2)
                if r.is_reverse == True:
                    direct1 = '-'
                else:
                    direct1 = '+'
                read_pos1 = r.query_alignment_end
                read_pos2 = s2
                DEL_len=(read_pos1-read_pos2)-(pos1-pos2)
                if DEL_len >= SV_LEN and direct2 == direct1:
                    long_DEL_record.append([chr1, min(pos1,pos2),max(pos1,pos2),DEL_len])

            elif chr1 == chr2:
                pos1 = r.reference_start
                pos2 = int(supple[1])
                direct2 = supple[2]
                cigar2 = supple[3]
                r2_len = 0
                i = 0
                while i < len(cigar2):
                    num = ''
                    while cigar2[i].isdigit():
                        num = num + cigar2[i]
                        i += 1
                    num = int(num)
                    if cigar2[i] == 'M' or cigar2[i] == 'D':
                        r2_len = r2_len + num
                    i += 1
                s2 = 0
                if cigar2[-1] == 'S' or cigar2[-1] == 'H':
                    s2 = num
                pos2 = pos2 + r2_len
                if r.is_reverse == True:
                    direct1 = '-'
                else:
                    direct1 = '+'
                read_pos2 = len(read) - s2
                read_pos1 = r.query_alignment_start
                DEL_len = (read_pos2 - read_pos1) - (pos2 - pos1)
                if DEL_len >= SV_LEN and direct2 == direct1:
                    long_DEL_record.append([chr1, min(pos1,pos2),max(pos1,pos2),DEL_len])

res=[]
if len(long_DEL_record):
    intervals = list(sorted(long_DEL_record))
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = intervals[0][3]
    read_num = 1
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][3]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr and \
                SV_len <= (high_rate * newSV_len) and SV_len >= (low_rate * newSV_len):
            low = intervals[i][1]
            high = intervals[i][2]
            read_num += 1
        else:
            res.append([chr, low, high, read_num,SV_len])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][3]
            read_num = 1
    res.append([chr, low, high, read_num,SV_len])
long_DEL=[]
for site in res:
    if site[3]>=3:
        long_DEL.append(site)
samfile.close()
#merge long DEL and short DEL.
DEL=short_DEL+long_DEL
res=[]
if len(DEL):
    intervals = list(sorted(DEL))
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = intervals[0][4]
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][4]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr and \
                SV_len <= (high_rate * newSV_len) and SV_len >= (low_rate * newSV_len):
            low = intervals[i][1]
            high = intervals[i][2]
        else:
            res.append([chr, low, high,int(SV_len)])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][4]
    res.append([chr, low, high,int(SV_len)])
samfile.close()

genes = Fasta(args.fasta)
candidate_res=[]
for site in res:
    if site[3]>=40 and site[3]<=200000 and (str(site[0]) in main_chr):
        chr=str(site[0])
        pos1 = int(site[1])
        pos2 = int(site[2])
        sv_len = int(site[3])
        ref = genes[chr][pos1:pos2]
        A_num = 0
        G_num = 0
        C_num = 0
        T_num = 0
        for base in ref:
            if base == 'a' or base == 'A':
                A_num += 1
            elif base == 'g' or base == 'G':
                G_num += 1
            elif base == 'c' or base == 'C':
                C_num += 1
            elif base == 't' or base == 'T':
                T_num += 1
        A_rate = float(A_num) / float(sv_len)
        G_rate = float(G_num) / float(sv_len)
        C_rate = float(C_num) / float(sv_len)
        T_rate = float(T_num) / float(sv_len)

        base_num = [['A', A_num], ['G', G_num], ['C', C_num], ['T', T_num]]
        base_num.sort(key=lambda x: x[1])
        ref_sign = [0] * len(ref)

        three_base = base_num[-1][0] * 10
        i = 0
        while i < (len(ref) - 10):
            if ref[i:i + 10] == three_base:
                ref_sign[i:i + 10] = '1111111111'
            i += 1
        base_sign = 0
        for site1 in ref_sign:
            if site1 == '1':
                base_sign += 1
        if len(ref) == 0:
            SR_rate == 0
        else:
            SR_rate = float(base_sign) / float(len(ref))

        if site[3]<=60:
            if (A_rate < 0.9 and G_rate < 0.9 and C_rate < 0.9 and T_rate < 0.9) and\
                    (SR_rate<0.4 and base_sign<20):
                candidate_res.append([str(site[0]),int(site[1]),int(site[2]),int(site[3]),0])
        else:
            candidate_res.append([str(site[0]),int(site[1]),int(site[2]),int(site[3]),0])

total_DEL=candidate_res

if len(total_DEL)!=0:
    DEL_final_res.append([str(total_DEL[0][0]), int(total_DEL[0][1]),int(total_DEL[0][2]),int(total_DEL[0][3]),'DEL'])
i=1
while i<len(total_DEL):
    if total_DEL[i-1][0]==total_DEL[i][0] and (int(total_DEL[i][1])-int(total_DEL[i-1][2]))<1000 and\
            int(total_DEL[i-1][3])<10000 and int(total_DEL[i][3])<10000:
        dif_rate = float(abs(total_DEL[i - 1][3] - total_DEL[i][3])) / float(min(total_DEL[i - 1][3], total_DEL[i][3]))
        if dif_rate < 0.2 and total_DEL[i - 1][3] > 200 and total_DEL[i][3] > 200:
            total_DEL[i - 1][4] = 2
            total_DEL[i][4] = 2
        else:
            total_DEL[i - 1][4] = 1
            total_DEL[i][4] = 1
            DEL_final_res.append(
                [str(total_DEL[i][0]), int(total_DEL[i][1]), int(total_DEL[i][2]), int(total_DEL[i][3]), 'DEL'])
    else:
        DEL_final_res.append(
            [str(total_DEL[i][0]), int(total_DEL[i][1]), int(total_DEL[i][2]), int(total_DEL[i][3]), 'DEL'])
    i+=1
#I add modified tep to modify the DEL pos. and add soft clip read step.



samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
INS_final_res=[]
#point= open('site2.txt', 'w')
record=[]
#detect the candidate DEL
for r in samfile.fetch():
  if r.mapq>20 and r.is_supplementary==False:
      chr = r.reference_name
      start_pos=r.reference_start
      cigar=r.cigarstring
      read_str=r.query_sequence
      i=0
      pos=start_pos
      pos1=0
      pos2=0
      read_i = 0
      while i<len(cigar):
          num=''
          while cigar[i].isdigit():
              num=num+cigar[i]
              i+=1
          num=int(num)
          if cigar[i]=='M':
              pos=pos+num
              read_i = read_i + num
          elif cigar[i]=='I':
              pos1 = pos
              pos2 = pos + num
              if num >= 30:
                  INS_seq = read_str[read_i:read_i + num]
                  record.append([chr, pos1, pos2, INS_seq])
              read_i = read_i + num
          elif cigar[i]=='D':
              pos=pos+num
          elif cigar[i]=='S':
              read_i = read_i+num
          i+=1
samfile.close()
#merge the candidate DEL
res=[]

if len(record):
    intervals = list(sorted(record))
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = len(intervals[0][3])
    SV_sequence = intervals[0][3]
    read_num = 1
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][2] - intervals[i][1]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr:
            low = min(intervals[i][1],low)
            high = max(intervals[i][2],high)
            SV_sequence = intervals[i][3]
            read_num += 1
        else:
            res.append([chr, low, high, read_num, SV_sequence])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = len(intervals[i][3])
            SV_sequence = intervals[i][3]
            read_num = 1
    res.append([chr, low, high, read_num, SV_sequence])
newres=[]
for site in res:
    if int(site[3])>=3 and len(site[4])>=30:
        newres.append([site[0],site[1],site[2],site[4]])
res_region=newres
#modified the DEL pos(For same read, there are two alignments.)

newDEL=[]
samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
for region in res_region:
    chr = region[0]
    sv_pos1 = region[1]
    sv_pos2 = region[2]
    read_del=[]
    zero_DEL_num = 0
    coverage=0
    for r in samfile.fetch(chr, sv_pos1, sv_pos2+1):
        read_tag=r.tags
        align_num=1
        for tag in r.tags:
            if str(tag[0])=="SA":
                sign="SA"
                supple_align=tag
                other_align = supple_align[1].split(";")
                other_align = other_align[:-1]
                align_num=align_num+len(other_align)

        if r.mapq > 20 and r.is_supplementary == False and \
                r.reference_start<(sv_pos1 - 1000) and r.reference_end>(sv_pos2 + 1000):
            coverage+=1
            range1=sv_pos1 - 1000
            range2 = sv_pos2 + 1000
            DEL_num=0
            DEL_len=0
            DEL_len_10bp = 0
            DEL_pos=[]
            start_pos = r.reference_start
            cigar = r.cigarstring
            read_str = r.query_sequence
            read_name=r.query_name
            i = 0
            pos = start_pos
            pos1 = 0
            pos2 = 0
            read_i = 0
            while i < len(cigar):
                num = ''
                while cigar[i].isdigit():
                    num = num + cigar[i]
                    i += 1
                num = int(num)
                if cigar[i] == 'M':
                    pos = pos + num
                    read_i = read_i + num
                elif cigar[i] == 'I':
                    pos1 = pos
                    pos2 = pos+num
                    INS_seq = read_str[read_i:read_i + num]
                    if pos1 > range1 and pos1 < range2:
                        if num >= 40:
                            DEL_pos.append([chr, pos1, pos1,num])
                            DEL_num += 1
                            DEL_len += num
                    read_i = read_i + num
                elif cigar[i] == 'D':
                    pos = pos + num
                elif cigar[i] == 'S':
                    read_i = read_i + num
                i += 1
            if DEL_num!=0:
                read_del.append([DEL_num, DEL_len, DEL_pos])
            else:
                zero_DEL_num+=1
            if coverage>150:
                break

    if len(read_del)==0 or float(len(read_del))/(float(len(read_del))+float(zero_DEL_num))<0.2:
        continue
    read_del = sorted(read_del, key=lambda x: x[1])
    noise_len = int(len(read_del) / 20)
    if noise_len != 0:
        read_del = read_del[noise_len:-noise_len]

    sv_len_num = []
    sv_len = read_del[0][1]
    read = read_del[0]
    num = 1
    i = 1
    while i < (len(read_del)):
        if read_del[i][1] == sv_len:
            if read_del[i][0] < read[0]:
                read = read_del[i]
            num += 1
        else:
            sv_len_num.append([sv_len, num, read])
            sv_len = read_del[i][1]
            read = read_del[i]
            num = 1
        i += 1
    sv_len_num.append([sv_len, num, read])

    cengci_cluster = []
    if len(sv_len_num) <= 2:
        cengci_cluster.append(sv_len_num)
    while len(sv_len_num) > 2:
        distance = []
        distance_rate = []
        i = 1
        while i < len(sv_len_num):
            node_dis = sv_len_num[i][0] - sv_len_num[i - 1][0]
            node_dis_rate = float(sv_len_num[i][0] - sv_len_num[i - 1][0]) / float(sv_len_num[i - 1][0])
            distance.append(node_dis)
            distance_rate.append(node_dis_rate)
            i += 1

        i = 1
        min_dis = distance[0]
        min_dis_rate = distance_rate[0]
        sign = 0
        while i < len(distance):
            if distance[i] < min_dis:
                min_dis = distance[i]
                min_dis_rate = distance_rate[i]
                sign = i
            i += 1

        if min_dis_rate<=0.2 or min_dis<=200:
            node_support_num = sv_len_num[sign][1] + sv_len_num[sign + 1][1]
            node_len = (float(sv_len_num[sign][0]) * sv_len_num[sign][1] + float(sv_len_num[sign + 1][0]) *
                        sv_len_num[sign + 1][1])
            node_len = node_len / node_support_num
            if sv_len_num[sign][2][0] < sv_len_num[sign + 1][2][0]:
                node_read = sv_len_num[sign][2]
            elif sv_len_num[sign][2][0] > sv_len_num[sign + 1][2][0]:
                node_read = sv_len_num[sign + 1][2]
            else:
                if abs(sv_len_num[sign][2][1] - node_len) < abs(sv_len_num[sign + 1][2][1] - node_len):
                    node_read = sv_len_num[sign][2]
                else:
                    node_read = sv_len_num[sign + 1][2]
            del sv_len_num[sign:sign + 2]
            sv_len_num.insert(sign, [node_len, node_support_num, node_read])
            sv_len_num_2 = copy.deepcopy(sv_len_num)
            cengci_cluster.append(sv_len_num_2)
        else:
            break

    if len(sv_len_num) == 1:
        chr = sv_len_num[0][2][2][0][0]
        del_pos1 = sv_len_num[0][2][2][0][1]
        del_pos2 = sv_len_num[0][2][2][-1][2]
        del_len = sv_len_num[0][0]
        newDEL.append([chr, del_pos1, del_pos2, sv_len_num[0][1], del_len])
    elif len(sv_len_num) == 2:
        cengci_cluster.pop()

        sv_len_num = sorted(sv_len_num, key=lambda x: x[1])
        sv_len_num.reverse()
        sv_len_num = sv_len_num[0:2]

        total_num = sv_len_num[0][1] + sv_len_num[1][1]
        support_read1=sv_len_num[0][1]
        support_read2=sv_len_num[1][1]
        rate1 = float(sv_len_num[0][1]) / float(total_num)
        rate2 = float(sv_len_num[1][1]) / float(total_num)
        del1=[sv_len_num[0][2][2][0][0],sv_len_num[0][2][2][0][1],sv_len_num[0][2][2][-1][2],sv_len_num[0][1],sv_len_num[0][0]]
        del2 = [sv_len_num[1][2][2][0][0], sv_len_num[1][2][2][0][1], sv_len_num[1][2][2][-1][2],sv_len_num[1][1],sv_len_num[1][0]]
        average_len = sv_len_num[0][1] * sv_len_num[0][0] + sv_len_num[1][1] * sv_len_num[1][0]
        average_len = float(average_len) / float(total_num)
        len_dif = abs(sv_len_num[1][0] - sv_len_num[0][0])
        while rate1<0.2 or rate2<0.2:
            if len(cengci_cluster)==0:
                rate1 = 0.5
                rate2 = 0.5
                len_dif = 0
                break
            cluster=sorted(cengci_cluster[-1], key=lambda x: x[1])
            total_num=cluster[-1][1]+cluster[-2][1]
            support_read1 = cluster[-1][1]
            support_read2 = cluster[-2][1]
            rate1=float(cluster[-1][1])/float(total_num)
            rate2 = float(cluster[-2][1]) / float(total_num)
            del1 = [cluster[-1][2][2][0][0], cluster[-1][2][2][0][1], cluster[-1][2][2][-1][2],cluster[-1][1],
                    cluster[-1][0]]
            del2 = [cluster[-2][2][2][0][0], cluster[-2][2][2][0][1], cluster[-2][2][2][-1][2], cluster[-2][1],
                    cluster[-2][0]]
            average_len = cluster[-1][1] * cluster[-1][0] + cluster[-1][1] * cluster[-1][0]
            average_len = float(average_len) / float(total_num)
            len_dif = abs(cluster[-1][0] - cluster[-2][0])
            cengci_cluster.pop()
            if cengci_cluster==0:
                break

        len_dif_value = 20 + float(0.005 * average_len)
        if len_dif > len_dif_value:
            if rate1 >= 0.2 and  support_read1 >= 3:
                newDEL.append(del1)
            if rate2 >= 0.2 and  support_read2 >= 3:
                newDEL.append(del2)
        else:
            if (rate1+rate2)>=0.2 and (support_read1+support_read2)>=3:
                newDEL.append(del1)

samfile.close()
short_DEL=newDEL

mini_len=300
#detect long DEL
samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
long_DEL_record=[]
for r in samfile.fetch():
    if r.mapq > 20 and r.is_supplementary == False:
        read_name=r.query_name
        sign = ""
        for tag in r.tags:
            if str(tag[0])=="SA":
                sign="SA"
                supple_align=tag
        # judge the two alignment
        if sign == "SA":
            # check primary alignment
            chr1 = r.reference_name
            cigar1 = r.cigarstring
            read = r.query_sequence
            other_align = supple_align[1].split(";")
            supple = other_align[0].split(",")
            i = 1
            while i < (len(other_align) - 1):
                m = other_align[i].split(",")
                if int(m[4]) > int(supple[4]):
                    supple = m
                i += 1
            chr2 = supple[0]
            # judge the primary alignment and supplementary alignment pos
            if r.reference_start <= int(supple[1]) and chr1 == chr2:
                pos1 = r.reference_end
                pos2 = int(supple[1])
                direct2 = supple[2]
                cigar2 = supple[3]
                s2 = cigar2.split("S")[0]
                if 'M' in s2:
                    s2 = 0
                else:
                    s2 = int(s2)
                if r.is_reverse == True:
                    direct1 = '-'
                else:
                    direct1 = '+'
                read_pos1 = r.query_alignment_end
                read_pos2 = s2
                INS_len = read_pos2 - read_pos1 - (pos2 - pos1)
                INS_seq = read[read_pos1:read_pos1+INS_len]
                if INS_len >= mini_len and direct2 == direct1:
                    long_DEL_record.append([chr1,min(pos1,pos2),max(pos1,pos2),INS_len])
            elif chr1 == chr2:
                pos1 = r.reference_start
                pos2 = int(supple[1])
                direct2 = supple[2]
                cigar2 = supple[3]
                r2_len = 0
                i = 0
                while i < len(cigar2):
                    num = ''
                    while cigar2[i].isdigit():
                        num = num + cigar2[i]
                        i += 1
                    num = int(num)
                    if cigar2[i] == 'M' or cigar2[i] == 'D':
                        r2_len = r2_len + num
                    i += 1
                s2 = 0
                if cigar2[-1] == 'S' or cigar2[-1] == 'H':
                    s2 = num
                pos2 = pos2 + r2_len
                if r.is_reverse == True:
                    direct1 = '-'
                else:
                    direct1 = '+'
                read_pos2 = len(read) - s2
                read_pos1 = r.query_alignment_start
                INS_len = read_pos1 - read_pos2-(pos1-pos2)
                INS_seq = read[read_pos2:read_pos2+INS_len]
                if INS_len >= mini_len and direct2 == direct1:
                    long_DEL_record.append([chr1,min(pos1,pos2),max(pos1,pos2),INS_len])

res=[]

if len(long_DEL_record):
    intervals = list(sorted(long_DEL_record))
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = intervals[0][3]
    #SV_sequence = intervals[0][4]
    read_num = 1
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][3]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr and \
                SV_len <= (high_rate * newSV_len) and SV_len >= (low_rate * newSV_len):
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][3]
            #SV_sequence = intervals[i][4]
            read_num += 1
        else:
            res.append([chr, low, high, read_num, SV_len])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][3]
            #SV_sequence = intervals[i][4]
            read_num = 1
    res.append([chr, low, high, read_num, SV_len])
long_DEL=[]
for site in res:
    if site[3]>=3:
        long_DEL.append(site)
samfile.close()
#merge long DEL and short DEL.
DEL=short_DEL+long_DEL

res=[]
if len(DEL):
    intervals = list(sorted(DEL))
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = intervals[0][4]
    #SV_sequence = intervals[0][5]
    read_num = 1
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][4]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr and \
                SV_len <= (high_rate * newSV_len) and SV_len >= (low_rate * newSV_len):
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][4]
            #SV_sequence = intervals[i][5]
            read_num += 1
        else:
            res.append([chr, low, high, SV_len])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][4]
            #SV_sequence = intervals[i][5]
            read_num = 1
    res.append([chr, low, high, SV_len])
samfile.close()

candidate_res=[]
for site in res:
    if site[3]>=40:
        candidate_res.append([str(site[0]),int(site[1]),int(site[2]),int(site[3]), 0])

total_DEL=candidate_res

if len(total_DEL)!=0:
    INS_final_res.append([str(total_DEL[0][0]),int(total_DEL[0][1]), int(total_DEL[0][2]), int(total_DEL[0][3]), 'INS'])
i=1
while i<len(total_DEL):
    if total_DEL[i-1][0]==total_DEL[i][0] and (int(total_DEL[i][1])-int(total_DEL[i-1][2]))<1000 and\
            int(total_DEL[i-1][3])<10000 and int(total_DEL[i][3])<10000:
        dif_rate = float(abs(total_DEL[i - 1][3] - total_DEL[i][3])) / float(min(total_DEL[i - 1][3], total_DEL[i][3]))
        if dif_rate < 0.2 and total_DEL[i - 1][3] > 200 and total_DEL[i][3] > 200:
            total_DEL[i - 1][4] = 2
            total_DEL[i][4] = 2
        else:
            total_DEL[i - 1][4] = 1
            total_DEL[i][4] = 1
            INS_final_res.append(
                [str(total_DEL[i][0]), int(total_DEL[i][1]), int(total_DEL[i][2]), int(total_DEL[i][3]), 'INS'])
    else:
        INS_final_res.append(
            [str(total_DEL[i][0]), int(total_DEL[i][1]), int(total_DEL[i][2]), int(total_DEL[i][3]), 'INS'])
    i+=1

# next..................................................................................

samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
INV_final_res=[]

SV_LEN=30
#detect long DEL
samfile=pysam.AlignmentFile(args.hg19_bam, 'rb')
long_DEL_record=[]
for r in samfile:
    if r.mapq > 20 and r.is_supplementary == False:
        sign = ""
        for tag in r.tags:
            if str(tag[0])=="SA":
                sign="SA"
                supple_align=tag
        # judge the two alignment
        if sign == "SA":
            # check primary alignment
            chr1 = r.reference_name
            cigar1 = r.cigarstring
            read = r.query_sequence
            other_align = supple_align[1].split(";")
            if r.is_reverse == True:
                direct1='-'
            else:
                direct1 = '+'

            supple=None
            MAPQ_value=0
            i=0
            while i < (len(other_align) - 1):
                m = other_align[i].split(",")
                if m[2]!= direct1 and int(m[4])>MAPQ_value:
                    supple = m
                    MAPQ_value=int(m[4])
                i+=1
            if supple==None:
                continue
            chr2 = supple[0]
            # judge the primary alignment and supplementary alignment pos
            if r.reference_start <= int(supple[1]) and chr1 == chr2:
                pos1 = r.reference_end
                pos2 = int(supple[1])
                cigar2 = supple[3]
                i=0
                middle_len=0
                while i < len(cigar2):
                    num = ''
                    while cigar2[i].isdigit():
                        num = num + cigar2[i]
                        i += 1
                    num = int(num)
                    if cigar2[i] == 'M' or cigar2[i] == 'D':
                        middle_len+=num
                    i+=1
                pos2 = pos2+middle_len

                DEL_len=pos2-pos1
                if DEL_len >= SV_LEN:
                    long_DEL_record.append([chr1, pos1, pos2])
            elif chr1 == chr2:
                pos1 = r.reference_start
                pos2 = int(supple[1])
                cigar2 = supple[3]
                DEL_len = pos1-pos2
                if DEL_len >= SV_LEN:
                    long_DEL_record.append([chr1, pos2, pos1])

res=[]
if len(long_DEL_record):
    intervals = list(sorted(long_DEL_record))
    '''for p in intervals:
        point.write(str(p[0])+':'+str(p[1])+'-'+str(p[2])+'\n')'''
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = intervals[0][2] - intervals[0][1]
    read_num = 1
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][2] - intervals[i][1]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr and \
                SV_len <= (1.2 * newSV_len) and SV_len >= (0.8 * newSV_len):
            low = intervals[i][1]
            high = intervals[i][2]
            read_num += 1
        else:
            res.append([chr, low, high, read_num])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][2] - intervals[i][1]
            read_num = 1
    res.append([chr, low, high, read_num])
long_DEL=[]
for site in res:
    if site[3]>=3:
        long_DEL.append(site)
samfile.close()
#merge long DEL and short DEL.
DEL=long_DEL

res=[]
if len(DEL):
    intervals = list(sorted(DEL))
    '''for p in intervals:
        point.write(str(p[0])+':'+str(p[1])+'-'+str(p[2])+'\n')'''
    chr = intervals[0][0]
    low = intervals[0][1]
    high = intervals[0][2]
    SV_len = intervals[0][2] - intervals[0][1]
    read_num = intervals[0][3]
    for i in range(1, len(intervals)):
        newSV_len = intervals[i][2] - intervals[i][1]
        if (high + merge_distance) >= intervals[i][1] and intervals[i][0] == chr and \
                SV_len <= (1.2 * newSV_len) and SV_len >= (0.8 * newSV_len):
            low = intervals[i][1]
            high = intervals[i][2]
            read_num += intervals[i][3]
        else:
            res.append([chr, low, high, read_num])
            chr = intervals[i][0]
            low = intervals[i][1]
            high = intervals[i][2]
            SV_len = intervals[i][2] - intervals[i][1]
            read_num = intervals[i][3]
    res.append([chr, low, high, read_num])

for site in res:
    if site[3]>20 and (site[2]-site[1])<100000:
       INV_final_res.append([str(site[0]), int(site[1]),int(site[2]),int(site[2]-site[1]),'INV'])
samfile.close()

localtime = time.asctime(time.localtime(time.time()))
print ('end',localtime)

localtime = time.strftime('%Y-%m-%d',time.localtime(time.time()))
w= open('result.txt', 'w')
w.write('##fileformat=VCFv4.2'+'\n')
w.write('##source=SVsearcher'+'\n')
w.write('##fileDate='+str(localtime)+'\n')

w.write('##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">'+'\n')
w.write('##ALT=<ID=DEL,Description="Deletion relative to the reference">'+'\n')
w.write('##ALT=<ID=INV,Description="Inversion of reference sequence">'+'\n')
w.write('##ALT=<ID=BND,Description="Breakend of translocation">'+'\n')

for chrom in genes.keys():
    w.write('##contig=<ID='+str(chrom) +',length='+str(len(genes[chrom]))+'>'+'\n')

w.write('##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variant">'+'\n')
w.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variant">'+'\n')
w.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'+'\n')
w.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'+'\n')
w.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">'+'\n')
w.write('##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">'+'\n')
w.write('##INFO=<ID=CILEN,Number=2,Type=Integer,Description="Confidence interval around inserted/deleted material between breakends">'+'\n')
w.write('##INFO=<ID=RE,Number=1,Type=Integer,Description="Number of read support this record">'+'\n')
w.write('##INFO=<ID=STRAND,Number=A,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">'+'\n')
w.write('##INFO=<ID=RNAMES,Number=.,Type=String,Description="Supporting read names of SVs (comma separated)">'+'\n')
w.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency.">'+'\n')
w.write('##FILTER=<ID=q5,Description="Quality below 5">'+'\n')
w.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'+'\n')

w.write('##CommandLine=SVsearcher '+str(args.hg19_bam)+' '+str(args.fasta)+'\n')
w.write('#CHROM	POS ID REF ALT QUAL FILTER INFO FORMAT NULL'+'\n')





final_sv_result=DEL_final_res+INS_final_res+INV_final_res
final_sv_result2=list(sorted(final_sv_result))
del_IDnum=0
ins_IDnum=0
inv_IDnum=0
for site in final_sv_result2:
    if str(site[4])=='DEL':
        del_IDnum+=1
        ID='SVsearcher.del.'+str(del_IDnum)
    elif str(site[4])=='INS':
        ins_IDnum+=1
        ID = 'SVsearcher.ins.' + str(ins_IDnum)
    elif str(site[4])=='INV':
        inv_IDnum+=1
        ID = 'SVsearcher.inv.' + str(inv_IDnum)

    w.write(str(site[0]) + ' ' + str(site[1]) + ' '+ ID+' '+str(genes[str(site[0])][site[1]:site[2]+1])+' . 60 PASS '+\
        'SVTYPE='+str(site[4])+'; '+'SVLEN='+str(site[3])+'; '+'END='+str(site[2])+'; GT:DR:DV:PL:GQ	./.:.:12:.,.,.:.'+'\n')
#I add modified tep to modify the DEL pos. and add soft clip read step.

# Домаћи задатак 2

#### Инсталирање пакета


```python
!pip install pysam
```

    Requirement already satisfied: pysam in /opt/conda/lib/python3.6/site-packages (0.10.0)



```python
#Importovanje AlignmentFile, AlignedSegment iz paketa
from pysam import AlignmentFile, AlignedSegment

```

#### Отварање САМ/БАМ датотеке


```python
samfile = AlignmentFile("../project-files/merged-tumor.bam", "rb")
```

#### Учитавање првог елемента


```python
s = samfile.head(1)
firstSegment: AlignedSegment = next(s)
print(firstSegment)
```

    C0HVYACXX120402:7:1207:5722:57044	1187	20	9483248	27	76M	20	9483381	76	TTTTCAAACAGTATCTATGCCTGCCAAATGTGAACATATAAAAAAAAACCAGAATGTGCCATTCTGATTTAAACTG	array('B', [28, 28, 27, 29, 31, 30, 31, 31, 29, 31, 35, 30, 29, 31, 34, 30, 29, 23, 41, 32, 20, 30, 29, 34, 34, 29, 30, 31, 30, 30, 30, 33, 33, 26, 39, 12, 25, 19, 32, 30, 35, 28, 35, 33, 23, 33, 35, 36, 30, 38, 33, 41, 34, 35, 31, 33, 23, 30, 30, 36, 27, 32, 29, 34, 35, 41, 33, 31, 33, 29, 32, 32, 31, 31, 31, 34])	[('XA', 'GL000217.1,-110754,76M,1;'), ('BD', 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'), ('MD', '76'), ('RG', '1'), ('BI', 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'), ('NM', 0), ('MQ', 27), ('AS', 76), ('XS', 71)]


#### Приказ поља објекта класе AligmentSegment


```python
print(dir(firstSegment))
#import inspect
#inspect.getmembers(firstSegment)
```

    ['__class__', '__copy__', '__deepcopy__', '__delattr__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__', '__pyx_vtable__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', 'aend', 'alen', 'aligned_pairs', 'bin', 'blocks', 'cigar', 'cigarstring', 'cigartuples', 'compare', 'flag', 'get_aligned_pairs', 'get_blocks', 'get_cigar_stats', 'get_overlap', 'get_reference_positions', 'get_reference_sequence', 'get_tag', 'get_tags', 'has_tag', 'infer_query_length', 'inferred_length', 'is_duplicate', 'is_paired', 'is_proper_pair', 'is_qcfail', 'is_read1', 'is_read2', 'is_reverse', 'is_secondary', 'is_supplementary', 'is_unmapped', 'isize', 'mapping_quality', 'mapq', 'mate_is_reverse', 'mate_is_unmapped', 'mpos', 'mrnm', 'next_reference_id', 'next_reference_name', 'next_reference_start', 'opt', 'overlap', 'pnext', 'pos', 'positions', 'qend', 'qlen', 'qname', 'qqual', 'qstart', 'qual', 'query', 'query_alignment_end', 'query_alignment_length', 'query_alignment_qualities', 'query_alignment_sequence', 'query_alignment_start', 'query_length', 'query_name', 'query_qualities', 'query_sequence', 'reference_end', 'reference_id', 'reference_length', 'reference_name', 'reference_start', 'rlen', 'rname', 'rnext', 'seq', 'setTag', 'set_tag', 'set_tags', 'tags', 'template_length', 'tid', 'tlen', 'tostring']


#### Провера флегова:


```python
for segment in samfile.head(5):
    print("----------------------")
    print("Flagovi: " + "0b{:012b} 0x{:03x}".format(segment.flag, segment.flag))
    sf = segment.flag

    if (sf & 0x1 > 0):
        print("read paired (0x1)")
    if (sf & 0x10 > 0):
        print("read reverse strand")
    if (sf & 0x40 > 0):
        print("first in pair")
    if (sf & 0x80 > 0):
        print("second in pair")
    if (sf & 0x400 > 0):
        print("read is PCR or optical duplicate (0x400)")
        
print("----------------------") 
```

    ----------------------
    Flagovi: 0b010010100011 0x4a3
    read paired (0x1)
    second in pair
    read is PCR or optical duplicate (0x400)
    ----------------------
    Flagovi: 0b000010100011 0x0a3
    read paired (0x1)
    second in pair
    ----------------------
    Flagovi: 0b000001100011 0x063
    read paired (0x1)
    first in pair
    ----------------------
    Flagovi: 0b000001100011 0x063
    read paired (0x1)
    first in pair
    ----------------------
    Flagovi: 0b000001100011 0x063
    read paired (0x1)
    first in pair
    ----------------------


## 2. део - Израчунавање


```python
mappedCount = 0
unmappedCount = 0
totalReadCount = 0

for read in samfile.fetch():
    totalReadCount+=1
    if read.is_unmapped:
        #read.flag & 0x4 != 0
        unmappedCount+=1
    else:
        mappedCount+=1
#print(totalReadCount)        
#print(mappedCount)
#print(unmappedCount)
print(f'Укупан број ридова: {totalReadCount}')
print(f'Број мапираних ридова: {mappedCount}')
print(f'Број немапираних ридова: {unmappedCount}')
```

    Укупан број ридова: 2921629
    Број мапираних ридова: 2903864
    Број немапираних ридова: 17765



```python
sumMappingQuality  = 0
mappingQualityZeroCount = 0
zeroAndMapped = 0
#minQual = 10000
#maxQual = -1
for read in samfile.fetch():
    #minQual = min(minQual, read.mapping_quality) #test
    #maxQual = max(maxQual, read.mapping_quality) #test
    sumMappingQuality += read.mapping_quality
    if read.mapping_quality == 0:
        mappingQualityZeroCount+=1
    if read.mapping_quality == 0 and not read.is_unmapped:#test
        zeroAndMapped += 1                                #test

averageMappingQuality = sumMappingQuality / totalReadCount
averageMappingQualityWithoutZero = sumMappingQuality / (totalReadCount - mappingQualityZeroCount)


print(averageMappingQuality, averageMappingQualityWithoutZero, mappingQualityZeroCount)#, zeroAndMapped
print(" ")
print(" ")
print("Просечан квалитет секвенцирања: " + "{:.2f}".format(averageMappingQuality))
print(
    "Просечан квалитет секвенцирања без ридова са квалитетом једнаким нули: " + 
    "{:.2f}".format(averageMappingQualityWithoutZero)
)
print(f'Број ридова са квалитетом секвенцирања једнаким нули {mappingQualityZeroCount}')
```

    55.91379158681681 58.446975510921106 126628
     
     
    Просечан квалитет секвенцирања: 55.91
    Просечан квалитет секвенцирања без ридова са квалитетом једнаким нули: 58.45
    Број ридова са квалитетом секвенцирања једнаким нули 126628


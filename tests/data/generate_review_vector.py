import sys,struct,pathlib
sys.dont_write_bytecode=True
# Supply the path to straw/C++/tests (independent Zstandard fixture helpers).
sys.path.insert(0, sys.argv[1])
from test_v10 import compress,pack,string
N=70000000
v=string('review')+pack('I',0)+pack('I',1)+string('chrA')+pack('Q',N)+pack('I',1)+pack('IBBHI',1,0,1,0,0xffffffff)+pack('I',0)+pack('I',1)+string('VC')
d=bytearray(88)+v
mp=len(d);d+=b'H10M'+pack('5I',1,0,0,1,0)+bytes(76)
bp=len(d)
raw=pack('5B3x4IQ2I',1,0,2,0,0,0,0,1,1,1,1,1)+b'\x00\x05'
b=b'H10B'+pack('BBHII',1,1,0,len(raw),0)+compress(raw);d+=b
ip=len(d);idx=b'H10I'+pack('IQII',2,40,1,0)+pack('IIQ',0,len(b),bp);d+=idx
d[mp+24:mp+100]=pack('4B3IB3xQQII2IQQ2I',0,0,1,0,0,1,0xffffffff,1,5,1,0x7fc00000,0x7fc00000,1000000,70,ip,40,1,0)
chunks=[]
for start in range(0,N,65536):
 n=min(65536,N-start);raw=pack('f',2.0)*n
 b=b'H10V'+pack('BBHII',1,0,0,len(raw),n)+compress(raw)
 chunks.append(pack('QIBBHQII',start,n,0,1,0,len(d),len(b),len(raw)));d+=b
ni=len(d);entry=pack('3IB3xIIQII',40+32*len(chunks),0,0,0,0,1,N,65536,len(chunks))+b''.join(chunks)
nvi=b'NVI0'+pack('III',1,1,0)+entry;d+=nvi
fp=len(d);footer=b'H10F'+pack('IQII',1,48,1,0)+pack('IIQQ',0,0,mp,100);d+=footer
d[:88]=b'HIC\0'+pack('I9QII',10,88+len(v),fp,48,ni,len(nvi),0,0,0,0,0,0)
p=pathlib.Path(__file__).with_name('large-vector-v10.hic');p.write_bytes(d);print('file bytes',len(d),'vector values',N)

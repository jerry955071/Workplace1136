from multiprocessing import Pool, TimeoutError, set_start_method, get_start_method
import time
import os
import random
import sys
import datetime as dt
import glob
import json

# class ComplexEncoder(json.JSONEncoder):
#     def default(self, obj):
#         typ = type(obj)
#         if typ in ['int', 'list', 'str', 'dict', 'float']:
#             return json.JSONEncoder.default(self, obj)
#         # if isinstance(obj, np.int64):
#         #     return int(obj)
#         # if isinstance(obj, KorkDay):
#         #     return str(obj)
#         return obj.__dict__


# class VO(object):
#     def __init__(self, code=0, status='OK'):
#         self.code = code
#         self.status = status
#         # self.payload = None
#         # self.data = None
#         self.misc = {}

#     def set_misc(self, obj):
#         self.misc = obj
#         return self

#     def set_err(self, so):
#         self.misc = so
#         return self

#     def tojson(self, indent=None):
#         return VO.toJson(self, indent=indent)

#     @staticmethod
#     def toJson(obj, indent=None):
#         try:
#             return json.dumps(obj, cls=ComplexEncoder, indent=indent, ensure_ascii=False)
#         except Exception as ex:
#             print(ex)
#             raise ex

#     @staticmethod
#     def fromJson(si):
#         try:
#             return json.loads(si)
#         except Exception as ex:
#             print(ex)
#             return {}

#     @staticmethod
#     def inFile(fname):
#         try:
#             with open(fname, 'r', encoding='utf-8') as fin:
#                 s = fin.read()
#             return s
#         except Exception as ex:
#             # s = traceback.format_exc()
#             # print(s)
#             print(ex)
#             return None


'''
for f in filename*.txt
do
    split -d -a1 -l10000 --additional-suffix=.txt "$f" "${f%.txt}-"
done
'''

def ParaRunner(pattern=r'*.fastq', outdir=None, forks=-1, **kwds):
    print('runner.enter.', forks)

    def deco(func):
        print('deco.enter', forks, pattern)

        def inner(*args, **kwds):
            # This is a MUST! Since we cannot modify the outter variables here.
            iforks = forks
            if iforks == -1:
                iforks = os.cpu_count()
            print('runner.inner.enter, forks=%s, pattern=%s' %
                  (iforks, pattern))
            if func is None:
                print('runner.inner.leave.NONE')
                return
            print('glob patterh=', pattern)
            lst = []
            for f in glob.glob(pattern):
                dic = {'infile': f, 'outfile': None}
                lst.append((f, None))
                # lst.append(dic)
                # func(dic)

            with Pool(processes=iforks) as pool:
                for r in pool.imap(func, lst):
                    print(r)
                print('fi')

            print('runner.inner.leave')

        return inner

    print('runner.leave')
    return deco

def ParaRunner2(func, pattern=r'*.fastq', outdir=None, forks=-1, **kwds):
    print('runner.enter.', forks)

    # Todo
    # if outdir is not None -> mkdir(outdir)

    # # This is a MUST! Since we cannot modify the outter variables here.
    # iforks = forks
    if forks == -1:
        forks = os.cpu_count()
    # print('runner.inner.enter, forks=%s, pattern=%s' %
    #         (iforks, pattern))
    if func is None:
        print('runner.inner.leave.NONE')
        return
    print('glob patterh=', pattern)
    lst = []
    for f in glob.glob(pattern):
        # fo = kutls.change_suffix(f, '.out')
        # if outdir is not None -> ofname = join(outdir, substitution(f))
        dic = {'infile': f, 'outfile': None}
        # lst.append((f, None))
        lst.append(dic)
        func(dic)

    with Pool(processes=forks) as pool:
        for r in pool.imap(func, lst):
            print(r)
        print('fi')

    print('runner.inner.leave')

def func(infile=None, outfile=None, *args, **kwds):
    if infile == None:
        print('Null input')
        return None, None

    if outfile == None:
        fout = sys.stdout
    else:
        try:
            fout = open(outfile, 'w')
        except Exception as ex:
            print(ex)
            return None, None
    try:
        with open(infile, 'r', encoding='UTF-8') as fin:
            for line in fin:
                info = line.strip()
                if info == '':
                    continue
                seq = fin.readline().strip()
                fin.readline().strip()
                qscore = fin.readline().strip()
                print(info[::-1], file=fout)
                print(seq[::-1], file=fout)
                print(qscore[::-1], file=fout)

    except Exception as ex:
        print(ex)
        return None, None

    finally:
        if fout is not sys.stdout:
            fout.close()

def proc2(tp):
    infile, outfile = tp
    # outfile = tp[1]
    print('process for (in,out)=(%s,%s).enter' % (infile, outfile))


# @ParaRunner(pattern=r'.in/*.fastq', forks=os.cpu_count())
def proc1(dic):
    outfile = dic['outfile']
    infile = dic['infile']
    print('process for (in,out)=(%s,%s).enter' % (infile, outfile))
    if infile == None:
        print('Null input')
        return None, None

    if outfile == None:
        fout = sys.stdout
    else:
        try:
            fout = open(outfile, 'w')
        except Exception as ex:
            print(ex)
            return None, None
    try:
        with open(infile, 'r', encoding='UTF-8') as fin:
            for line in fin:
                # seq = fin.readline()
                # seq = fin.readline()
                # seq = fin.readline()
                # seq = fin.readline()
                # time.sleep(random.random()*1.0)
                info = line.strip()
                if info == '':
                    continue
                if info[0:1] != '@':
                    print('you have a trouble!')
                # continue
                seq = fin.readline().strip()
                fin.readline().strip()
                qscore = fin.readline().strip()
                # dx = { 'info': info, 'seq': seq, 'qscore': qscore}
                # s = VO.toJson(dx)
                # print(s)
                # dy = VO.fromJson(s)
                print(info[::-1], file=fout)
                print(seq[::-1], file=fout)
                print(qscore[::-1], file=fout)

    except Exception as ex:
        print(ex)
        return None, None

    finally:
        if fout is not sys.stdout:
            fout.close()
        print('process for (in,out)=(%s,%s).leave' % (infile, outfile))
        return infile+'.fin'

# func('.in/001.fastq')


def f(x):
    t = 0.1 + 3*random.random()
    pid = os.getpid()
    print('#%d.input=%s,to.sleep=%f' % (pid, x, t))
    # print('mypid=%d, parent=%d' % (os.getpid(), os.getppid()))
    time.sleep(t)
    # ret
    print('#%d,to.return=%f,for.input=%s' % (pid, t, x))
    return t
    # print('input is', x, ', do output', x*x)
    # return x*x


def decoB(name: str = None, description: str = None, times: int = 2):
    print('decoB.%s enter' % name)

    def _pipeline(func, *args, **kwds):
        print('decoB.inner.%s enter' % name)
        func(*args, **kwds)
        # print('split %s into %d files' % (ifname, times))
        # for i in range(times):
        #     print('decoB.inner.%s callout.%d' % (name, i))
        #     func()
        print('decoB.inner.%s leave' % name)
        return
    print('decoB.%s leaves' % name)
    return _pipeline


# print('---a')
# @decoB
# def main2(*args, **kwds):
#     print('main2.called w/', args, kwds)

# main()
# main2()
# def get_files(tgtdir, pattern, reverse=False):
#         lst = sorted([str(e.resolve()) for e in Path(
#             tgtdir).glob(pattern)], reverse=reverse)
#         return lst

# import functools
def Runner(pattern=r'*.fastq', outdir=None, forks=9, **kwds):
    print('runner.enter.', forks)

    def deco(func):
        # forks = for
        print('deco.enter', forks, pattern)

        # @functools.wraps(inner)
        def inner(*args, **kwds):
            # forks = forks
            # if forks == 0:
            #     forks = os.cpu_count()
            print('runner.inner.enter, forks=', forks)
            dd = forks
            if dd == 0:
                dd = os.cpu_count()
            if dd > os.cpu_count():
                dd = os.cpu_count()
            print('forks=', dd)
            if func is None:
                print('runner.inner.leave.NONE')
                return
            print('glob patterh=', pattern)
            for f in glob.glob(pattern):
                func(infile=f, outfile=None)
            print('runner.inner.leave')

        return inner

    print('runner.leave')
    return deco


# @Runner(forks=99, xforks=os.cpu_count(), pattern=r'.in\*.fastq')
# def mainx(*args, **kwds):
#     print('xmain.called w/', args, kwds)


# mainx()


# @Runner(pattern='zz*')


# def f2(prms):
#     proc1(prms)

# f2()
# Runner(func=func)
# Runner(r'.in\*.fastq', func=func)
# Runner(r'.in/*.fastq', func=func)
# Runner(r'.in/[\d]{3}.fast*')
# sys.exit(0)
if __name__ == '__main__':
    # fx = decoB('nano.flv')
    # fx(main)
    # main(ifname='Zoo')
    # main2(main, ifname='Koo')

    # ParaRunner(pattern=r'.in/00*.fastq', forks=3)(proc1)()
    ParaRunner2(proc1, pattern=r'.in/00*.fastq', forks=3)
    
    # proc1()

    # Runner(forks=9)(f2)

    # @Runner(forks=os.cpu_count(), pattern=r'.in\*.fastq')
    # def main(*args, **kwds):
    #     print('main.called w/', args, kwds)

    # i = get_start_method()
    # set_start_method('fork')
    print('the host cpu core is', os.cpu_count())
    with Pool(processes=os.cpu_count()) as pool:

        # print "[0, 1, 4,..., 81]"
        # print(pool.map(f, range(10)))

        # print same numbers in arbitrary order
        # for i in pool.imap_unordered(f, [(3, 3), {'key': 'waht'}, dt.datetime.now()]):
        #     print(i)

        for i in range(100):
            time.sleep(0.2)
            fname = 'nano'+str(i)
            pool.imap(f, [fname])
            # for n in pool.imap(f, [(i, i)]):
            #     print(n)

        for i in range(100):
            time.sleep(0.2)
            fname = 'nano'+str(i+10000)
            pool.imap(f, [fname])
            # pool.imap_unordered(f, [(i+9, i+9)])
            # for n in pool.imap_unordered(f, [(i+9, i+9)]):
            #     print(n)

        for n in pool.imap_unordered(f, [(-1, -1)]):
            print(n)
        print('fi')


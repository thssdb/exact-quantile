import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.longs.LongLongPair;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import it.unimi.dsi.util.XoRoShiRo128PlusRandom;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class MainForExactSeparation {
    static int mergeBufferRatio = 5;
    static boolean DEBUG_PRINT=false;
    int dataType=-233;
    static int startType=5,endType=5;
    static int pageN=4096, memTableN=262144;
    static int N = (int) (1e8+1e7)/*(1e9 + 1e7 - 262144)*/, pageNum = N / pageN; // CHECK IT
    long seqN,unSeqN;
    public static int TEST_CASE = 1; // CHECK IT
    static double[] a;
    static ObjectArrayList<SSTable> seqSST,unSeqSST;
     int RESULT_LINE = 0;
    public XoRoShiRo128PlusRandom random = new XoRoShiRo128PlusRandom(233);
    static long[] aa;
    static long[] bb;
    static IntArrayList cc;
    static double mu, sig;


    private static class PageStat {
        LongKLLSketch sketch;
        int count=0;
        long minT,maxT,minV,maxV;
//        boolean overlappedWithPage=false;
        boolean overlappedWithPoint=false;
        public void add(long t,long longV) {
            sketch.update(longV);
            minT = Math.min(minT, t);
            maxT = Math.max(maxT, t);
            minV = Math.min(minV, longV);
            maxV = Math.max(maxV, longV);
            count++;
        }
    }

    private static class Page {
        PageStat stat=new PageStat();
        ObjectArrayList<LongLongPair> data=new ObjectArrayList<>(pageN);
        public Page(int SummaryByte){
            stat.sketch = new LongKLLSketch(pageN, pageN*8, SummaryByte);
            stat.minT=stat.minV=Long.MAX_VALUE;
            stat.maxT=stat.maxV=Long.MIN_VALUE;
        }
        public void add(long t,long longV) {
            stat.add(t,longV);
            data.add(LongLongPair.of(t,longV));
        }
        public void end(){stat.sketch.compactBeforeSerialization();}
    }
    private static class SSTable {
        int PageSummaryByte;
        ObjectArrayList<Page> pageList = new ObjectArrayList<>();
        int count=0;
        long minT,maxT,minV,maxV;
        public SSTable(int SummaryByte){
            PageSummaryByte=SummaryByte;
            minT=minV=Long.MAX_VALUE;
            maxT=maxV=Long.MIN_VALUE;
        }
        public SSTable(int SummaryByte,ObjectArrayList<LongLongPair> dataList){
            PageSummaryByte=SummaryByte;
            minT=minV=Long.MAX_VALUE;
            maxT=maxV=Long.MIN_VALUE;
            build(dataList);
        }
        public void build(ObjectArrayList<LongLongPair> dataList) {
            dataList.sort(Comparator.comparingLong(LongLongPair::firstLong));
            Page cntPage = new Page(PageSummaryByte);
            for(LongLongPair data:dataList){
                cntPage.add(data.firstLong(),data.secondLong());
                if(cntPage.stat.count==pageN){
                    cntPage.end();
                    pageList.add(cntPage);
                    cntPage=new Page(PageSummaryByte);
                }
            }
            if(cntPage.stat.count>0){
                cntPage.end();
                pageList.add(cntPage);
            }
            for(Page page:pageList){
                count+=page.stat.count;
                minT=Math.min(minT,page.stat.minT);
                maxT=Math.max(maxT,page.stat.maxT);
                minV=Math.min(minV,page.stat.minV);
                maxV=Math.max(maxV,page.stat.maxV);
            }
        }
    }


    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = //longToResult(i&4095);
                    //Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
                    random.nextGaussian();
//                    i-0.5*Math.floor(i/pageN)*pageN;
//                    i;
            return;
        }
        if (dataType == 5) {
            for (int i = 0; i < N; i++)
                a[i] = random.nextGaussian();
            return;
        }
        BufferedReader reader = null;
        if (dataType == 1)
            reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
        if (dataType == 2)
            reader = new BufferedReader(new FileReader(new File("2_SpacecraftThruster.txt")));
        if (dataType == 3)
            reader = new BufferedReader(new FileReader(new File("3_taxipredition8M.txt")));
        if (dataType == 4)
            reader = new BufferedReader(new FileReader(new File("4_wh.csv")));
        assert reader != null;
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
    }

    public void prepareDisorder() {
        random = new XoRoShiRo128PlusRandom(233);
        aa = new long[N];
        bb = new long[N];
        cc = new IntArrayList(N);
        for (int i = 0; i < N; i++) {
            cc.add(i);
            aa[i] = i;
            bb[i] = Math.round(i + Math.exp(mu + sig * random.nextGaussian()));
        }
        cc.sort((x, y) -> (Long.compare(bb[x], bb[y])));
    }


    public void prepareWorker(int maxSeriesByte) {
        ObjectArrayList<LongLongPair> seqTVList = new ObjectArrayList<>(memTableN);
        ObjectArrayList<LongLongPair> unSeqTVList = new ObjectArrayList<>(memTableN);
        seqSST = new ObjectArrayList<>();
        unSeqSST=new ObjectArrayList<>();
        long lastFlushSeqEndTime = Long.MIN_VALUE;
        LongArrayList unSeqT = new LongArrayList();
        for(int i=0;i<N;i++){
            int index = cc.getInt(i);
            long cntT = aa[index];
            long cntV = dataToLong(a[index]);
            if(cntT<=lastFlushSeqEndTime){
                unSeqT.add(cntT);
                unSeqTVList.add(LongLongPair.of(cntT,cntV));
                if(unSeqTVList.size()==memTableN){
                    unSeqSST.add(new SSTable(maxSeriesByte,unSeqTVList));
                    unSeqTVList.clear();
                }
            }else {
                seqTVList.add(LongLongPair.of(cntT,cntV));
                if(seqTVList.size()==memTableN){
                    SSTable tmp = new SSTable(maxSeriesByte,seqTVList);
                    seqSST.add(tmp);
                    seqTVList.clear();
                    lastFlushSeqEndTime=tmp.pageList.get(tmp.pageList.size()-1).stat.maxT;
                }
            }
        }
        if(unSeqTVList.size()>0){ // then seqTVList.size()>0, since N=pageN*?
            unSeqSST.add(new SSTable(maxSeriesByte,unSeqTVList));
            unSeqTVList.clear();
        }
        if(seqTVList.size()>0){
            seqSST.add(new SSTable(maxSeriesByte,seqTVList));
            seqTVList.clear();
        }
        unSeqSST.sort(Comparator.comparingLong(x -> x.minT));

        unSeqT.sort(Long::compare);
        int cntUnSeqID=0;
        for(SSTable sst:seqSST)
            for(Page page:sst.pageList){
                while(cntUnSeqID<unSeqT.size()&&unSeqT.getLong(cntUnSeqID)<page.stat.minT)cntUnSeqID++;
                if(cntUnSeqID<unSeqT.size()&&unSeqT.getLong(cntUnSeqID)<=page.stat.maxT)page.stat.overlappedWithPoint=true;
                while(cntUnSeqID<unSeqT.size()&&unSeqT.getLong(cntUnSeqID)<=page.stat.maxT)cntUnSeqID++;
            }
        long seqN=0,unSeqN=0;
        for(SSTable sst:seqSST)seqN+=sst.count;
        for(SSTable sst:unSeqSST)unSeqN+=sst.count;
        System.out.println("\t\tseqSST num:\t"+seqSST.size()+"\t\tseqN:\t"+seqN+"\tseqNRatio:\t"+1.0*seqN/N);
        System.out.println("\t\tunSeqSST num:\t"+unSeqSST.size()+"\t\tunSeqN:\t"+unSeqN+"\tunSeqNRatio:\t"+1.0*unSeqN/N);
    }

    public int getValueActualRank(double[] sortedA, int queryN, double v) { // number of elements <= v
        int L = 0, R = queryN - 1;
        while (L < R) {
            int mid = (L + R + 1) >>> 1;
            if (v < sortedA[mid]) R = mid - 1;
            else L = mid;
        }
        return L;
    }


    private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }

    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }

    private boolean inInterval(long x, long y, long L, long R) {
        return x >= L && y <= R;
    }

    private boolean inInterval(long x, long L, long R) {
        return x >= L && x <= R;
    }
    private boolean overlap(long x, long y, long L, long R) { // [L,R]
        return !(y < L || x > R);
    }

    private boolean overlapInterval(long x, long y, long L, long R) { // [L,R]
        return !inInterval(x, y, L, R) && !(y < L || x > R);
    }

    private void mergeBuffered(ObjectArrayList<PageStat>buffer,KLLSketchLazyExact cntWorker){
        if(buffer.isEmpty())return;
        List<KLLSketchForQuantile> sketchList = new ArrayList<>();
        LongArrayList minVList=new LongArrayList(),maxVList=new LongArrayList();
        for(PageStat pageStat:buffer) {
            sketchList.add(pageStat.sketch);
            minVList.add(pageStat.minV);
            maxVList.add(pageStat.maxV);
        }
        cntWorker.mergeWithTempSpace(sketchList,minVList,maxVList);
        buffer.clear();
    }
    private void addPageStat(PageStat pageStatToAdd,ObjectArrayList<PageStat>buffer,int maxBufSize,KLLSketchLazyExact cntWorker){
        buffer.add(pageStatToAdd);
        if(buffer.size()>=maxBufSize){
            mergeBuffered(buffer,cntWorker);
        }
    }


    public void testError(int queryN, int maxMemoryByte,int pageSketchByte,boolean ignoreUpdate,boolean noUpdate) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
//        DecimalFormat fnum = new DecimalFormat("#0.00");
        long full_time = 0, merge_page_time = 0;
        double avg_iteration=0;
//        double[] query_a = new double[queryN];

        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        double ALLPrSum=0,ALLPrCount=0,FailCount=0,IterCount=0;
        double IOCost=0;

        for (int T = 0; T < TEST_CASE; T++) {
//            prepareWorker(pageSketchByte);
//            if(T%100==0) System.out.println("\tTEST_ID:"+T);
            int L = 0, R = queryN;
//            int L=LL[T],R=RR[T];

//            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
//            Arrays.sort(query_a);

            double q_tmp = 1e-2, q_start = q_tmp, q_end = 1 - q_tmp;
                q_start=0.5;q_end=0.5;
            double q_add = 1e-2, q_count = Math.floor((q_end - q_start + 1e-10) / q_add) + 1;
            double ratioPerQuery = 1.0 / (q_count * TEST_CASE);

            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
//                if(q<0.189||q>0.194)continue;
                int query_rank1 = (int) Math.floor(q * (queryN-1)+1),query_rank2 = (int) Math.ceil(q * (queryN-1)+1);
//                System.out.println("\n\t\t\t\t\t\tq rank1,2:"+q+"\t"+query_rank1+","+query_rank2);
                long last_n=queryN;
                double[] deterministic_result,iterate_result;
                iterate_result = new double[]{-Double.MAX_VALUE,Double.MAX_VALUE};
                deterministic_result = iterate_result;

                int MMP=0;
                while(deterministic_result[0]<deterministic_result[1]&&deterministic_result.length!=3) {
                    if(++MMP>15){
                        System.out.println("\t\t\t\t\t\titerate fail.\titer:"+MMP+"\t\tq:"+q);
                        return;
                    }
//                    System.out.println("\n\tnew iter.??\t"+iterate_result[0]+"..."+iterate_result[1]);
                    avg_iteration+=ratioPerQuery;
                    KLLSketchLazyExact cntWorker = new KLLSketchLazyExact(maxMemoryByte);
//                    FastKLLSketchLazyForBound cntWorker = new FastKLLSketchLazyForBound(maxMemoryByte);

                    IterCount+=1;
                    double valL=iterate_result[0],valR=iterate_result[1];
                    int CountOfLessThanValL=0;
                    int IgnoredUnSeqN=0;

//            System.out.println("\t\t\t"+posL+"\t"+posR);

                    if(deterministic_result[0]==-Double.MAX_VALUE&&deterministic_result[1]==Double.MAX_VALUE) {
                        cntWorker = new KLLSketchLazyExact(maxMemoryByte*(mergeBufferRatio-1)/mergeBufferRatio);
                        ObjectArrayList<PageStat> statBuffer=new ObjectArrayList<>();
                        int maxBufSize=maxMemoryByte/mergeBufferRatio/pageSketchByte;

                        for(SSTable sst:seqSST)
                            if(overlap(sst.minT,sst.maxT,L,R-1)) {
                                for (Page page : sst.pageList)
                                    if (inInterval(page.stat.minT, page.stat.maxT, L, R-1)&&(ignoreUpdate||!page.stat.overlappedWithPoint))
                                        addPageStat(page.stat, statBuffer, maxBufSize, cntWorker);
                                    else if (overlap(page.stat.minT, page.stat.maxT, L, R-1)) {
                                        IOCost+=page.stat.count;
                                        for(LongLongPair data:page.data)
                                            if(inInterval(data.firstLong(),L,R-1))
                                                cntWorker.update(data.secondLong());
                                    }
                            }
                        for(SSTable sst:unSeqSST)
                            if(overlap(sst.minT,sst.maxT,L,R-1)) {
                                for (Page page : sst.pageList)
                                    if (inInterval(page.stat.minT, page.stat.maxT, L, R-1)&&(ignoreUpdate)) {
                                        IgnoredUnSeqN+=page.stat.count;
                                        addPageStat(page.stat, statBuffer, maxBufSize, cntWorker);
                                    }
                                    else if (overlap(page.stat.minT, page.stat.maxT, L, R-1)) {
                                        IOCost+=page.stat.count;
                                        for(LongLongPair data:page.data)
                                            if(inInterval(data.firstLong(),L,R-1)) {
                                                IgnoredUnSeqN++;
                                                cntWorker.update(data.secondLong());
                                            }
                                    }
                            }
                        mergeBuffered(statBuffer,cntWorker);
                        if(ignoreUpdate&&!noUpdate)System.out.println("\t\t\ttemp_IgnoredUnSeqN:\t"+IgnoredUnSeqN);
                    }else {
//                        IOCost+=R-L;
//                        for (int i = L; i < R; i++) {
//                            if (a[i] >= valL && a[i] <= valR) cntWorker.update(dataToLong(a[i]));
//                            else if (a[i] < valL) CountOfLessThanValL++;
//                        }


                        for(SSTable sst:seqSST)
                            if(overlap(sst.minT,sst.maxT,L,R-1)) {
                                for (Page page : sst.pageList)
                                    if (overlap(page.stat.minT, page.stat.maxT, L, R-1)) {
                                        IOCost+=page.stat.count;
                                        for(LongLongPair data:page.data)
                                            if(inInterval(data.firstLong(),L,R-1)) {
                                                double douV=longToResult(data.secondLong());
                                                if (douV >= valL && douV <= valR) cntWorker.update(dataToLong(douV));
                                                else if (douV < valL) CountOfLessThanValL++;
                                            }
                                    }
                            }
                        for(SSTable sst:unSeqSST)
                            if(overlap(sst.minT,sst.maxT,L,R-1)) {
                                for (Page page : sst.pageList)
                                    if (overlap(page.stat.minT, page.stat.maxT, L, R-1)) {
                                        IOCost+=page.stat.count;
                                        for(LongLongPair data:page.data)
                                            if(inInterval(data.firstLong(),L,R-1)) {
                                                double douV=longToResult(data.secondLong());
                                                if (douV >= valL && douV <= valR) cntWorker.update(dataToLong(douV));
                                                else if (douV < valL) CountOfLessThanValL++;
                                            }
                                    }
                            }
                    }
//                    cntWorker.show();
//                    cntWorker.showNum();
//                    cntWorker.showCompact();
//                    if(MMP==1&&T==0&&q==q_start)cntWorker.showCompact();
//                    if(cntWorker.getN()==cntWorker.getNumLen())System.out.println("\t\t\tfinal iter\tspaceRate:"+1.0*cntWorker.getN()/(maxMemoryByte/8));
//                    cntWorker.show();
//                    System.out.println("\t\t\t\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntSketch_N:"+cntWorker.getN());
                    int cntRank1=query_rank1-CountOfLessThanValL;
                    int cntRank2=query_rank2-CountOfLessThanValL;

                    if(ignoreUpdate&&IgnoredUnSeqN>0&&!noUpdate){
                        cntRank1 = Math.max(1,cntRank1-IgnoredUnSeqN/2);
                        cntRank2 = Math.min((int)cntWorker.getN(),cntRank2+IgnoredUnSeqN/2);
                        System.out.println("\t\t\ttemp_cntRank:"+cntRank1+" "+cntRank2);
                    }

//                    System.out.println("\t\t\t\t\t\tcntRank:"+cntRank1+" "+cntRank2);
                    if(cntRank1<=0||cntRank2>cntWorker.getN()){ // iteration failed.
                        if(DEBUG_PRINT)System.out.println("\t\t\t\t\t\titerate fail."+"\t\tcntIter:"+MMP);
                        if(cntRank1<=0)
                            iterate_result = new double[]{deterministic_result[0],iterate_result[0]};
                        else iterate_result = new double[]{iterate_result[1],deterministic_result[1]};
                        deterministic_result = iterate_result;
                        if(deterministic_result[0]==deterministic_result[1])continue;
                        FailCount+=1;
//                        IterCount-=1;
                        continue;
                    }
                    if(DEBUG_PRINT)System.out.println("\t\t\t\t\t\titerate success."+"\t\tcntN:"+cntWorker.getN()+"\t\tcntIter:"+MMP);
                    deterministic_result = cntWorker.findResultRange(cntRank1,cntRank2,1.0);
                    last_n=cntWorker.getN();
                    if(deterministic_result.length==3){
                        iterate_result = deterministic_result;
                        break;
                    }

                    FindBestPrHelper findPr = new FindBestPrHelper(DEBUG_PRINT);
                    double bestPr=findPr.findBestPr(cntWorker,maxMemoryByte/8,cntRank1,cntRank2,deterministic_result);

                    ALLPrSum+=bestPr;ALLPrCount+=1;
                    iterate_result = cntWorker.findResultRange(cntRank1,cntRank2,bestPr);
//                    System.out.println("\t\t\t\tbestPr:"+bestPr+"\t\t\tcntN:"+cntWorker.getN());
//                    System.out.println("\t\t\t\t\tcntL,R:"+iterate_result[0]+","+iterate_result[1]+"\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntRank:"+cntRank1+" "+cntRank2);
                }
                double exact_quantile_v = (deterministic_result[0]+deterministic_result[1])*0.5;
                if(DEBUG_PRINT)System.out.println("FINISH CALC\t\tq:"+q+"\texact_quantile="+exact_quantile_v+"\t\titer:\t"+MMP+"\tIOCost:\t"+IOCost+"\n");
//                if(Math.abs(exact_quantile_v-N*q)>1.0){
//                    System.out.println("!!!!!!!!!!ERROR!!!!");
//                    return;
//                }
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("FindFixPR!\tTEST_CASE="+TEST_CASE+"\tDATASET:"+dataType+"\tqueryN:\t" + queryN+"\tmemory:\t" + maxMemoryByte
            +"\t\tmu:\t"+mu+"\tsig:\t"+sig+"\t\tignoreUpdate:\t"+ignoreUpdate+"\t\tnoUpdate:\t"+noUpdate
            +"\t\tavg_iteration:\t"+avg_iteration+"\t\tavgIOCost:\t"+IOCost/TEST_CASE+"\t\tavgPrChose:"+ALLPrSum/ALLPrCount+"\t\tavgFailRate:"+(FailCount/IterCount)+"\tusedNorDistri:\t"+KLLSketchLazyExact.DEBUG_COUNT_USENORMAL+"\tnot_cached_calc_pr_err:"+KLLSketchLazyExact.DEBUG_COUNT_USEPRERR);
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactSeparation main;

        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main = new MainForExactSeparation();
            main.prepareA(dataType);
            double mu = 3.0;
            for (double sig:new double[]{0,0.6,1.2,1.8,2.1,2.3,2.5,2.65,2.8,2.9,3.0}) {
                MainForExactSeparation.mu = mu;
                MainForExactSeparation.sig = sig;
                main.prepareDisorder();
                main.prepareWorker(256);;
                for (int queryN : new int[]{(int)5e7})
                    for (int query_mem : new int[]{1024*1024})
                        for (int page_seri : new int[]{256}) {
                            main.testError(queryN, query_mem, page_seri, true,true);
                            main.testError(queryN, query_mem, page_seri, true,false);
                            main.testError(queryN, query_mem, page_seri, false,false);
                        }
                System.out.println("\n---------------------------\n");
            }
//            for (double sig:new double[]{3.0}) {
//                MainForExactSeparation.mu = mu;
//                MainForExactSeparation.sig = sig;
//                main.prepareDisorder();
//                main.prepareWorker(256);;
//                for (int queryN : new int[]{(int)5e7})
//                    for (int query_mem : new int[]{1024*1024,1024*1024*2,1024*1024*4,1024*1024*6,1024*1024*8})
//                        for (int page_seri : new int[]{256}) {
//                            main.testError(queryN, query_mem, page_seri, true,true);
//                            main.testError(queryN, query_mem, page_seri, true,false);
//                            main.testError(queryN, query_mem, page_seri, false,false);
//                        }
//                System.out.println("\n---------------------------\n");
//            }
        }
        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}

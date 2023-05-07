import it.unimi.dsi.fastutil.doubles.DoubleArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.Random;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForExactLazyKLLFindPr {
    static boolean DEBUG_PRINT=false;
    int dataType=-233;
    static int startType=5,endType=5;
    static int pageN = 8192;
    static int N = pageN * 6712/*1000000*/, pageNum = N / pageN; // CHECK IT
    public static int TEST_CASE = 2048; // CHECK IT
    static double[] a;
    static DoubleArrayList prList=null;
    static DoubleArrayList simuPrList=null;
//    static double[] pageMinV,pageMaxV;
    static int compaction_level;
    ArrayList<String> time_result = new ArrayList<>();
    int RESULT_LINE = 0;
    Random random = new Random(233);

    static int DEBUG_COUNT_SIMULATION=0;

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = //longToResult(i&4095);
                    //Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
                    random.nextDouble();
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

//    public void prepareWorker(int maxSeriesByte) {
//        pageMinV = new double[pageNum];pageMaxV = new double[pageNum];
//
//        int enoughMemByte = pageN * 9;
//        for (int i = 0; i < pageNum; i++) {
//            pageMinV[i] = Long.MAX_VALUE;
//            pageMaxV[i] = Long.MIN_VALUE;
//            for (int j = 0; j < pageN; j++) {
//                double v =a[i*pageN+j];
//                pageMinV[i] = Math.min(pageMinV[i],v);
//                pageMaxV[i]=Math.max(pageMaxV[i],v);
//            }
////            worker.show();
//        }
//    }

    public int getValueActualRank(double[]sortedA,int queryN, double v){ // number of elements <= v
        int L=0,R=queryN-1;
        while(L<R){
            int mid=(L+R+1)>>>1;
            if(v<sortedA[mid])R=mid-1;
            else L=mid;
        }
        return L;
    }
    public int getValueLessThan(double[]sortedA,int queryN, double v){ // number of elements <= v
        int L=0,R=queryN-1;
        while(L<R){
            int mid=(L+R+1)>>>1;
            if(sortedA[mid]<v)L=mid;
            else R=mid-1;
        }
        return sortedA[L]<v?L:L-1;
    }
    public int getDeltaRank(double[]sortedA,int queryN, double v,int targetRank){
        int rank_L = getValueLessThan(sortedA,queryN,v)+1;
        int rank_R = getValueActualRank(sortedA,queryN,v);
//        System.out.println("\t\t\t"+targetRank+"\t\tresultLR:"+rank_L+"..."+rank_R+"\t\tresV:"+v);
        if(targetRank>=rank_L&&targetRank<=rank_R)return 0;
        else return targetRank<rank_L?(targetRank-rank_L):(targetRank-rank_R);
    }
    private long dataToLong(double data) {
        long result = Double.doubleToLongBits((double) data);
        return data >= 0d ? result : result ^ Long.MAX_VALUE;
    }
    private double longToResult(long result) {
        result = (result >>> 63) == 0 ? result : result ^ Long.MAX_VALUE;
        return Double.longBitsToDouble(result);
    }

    private double estimateIteration(int n, int maxMemoryNum){
        if(n<=0)return 0;
        if(n<=maxMemoryNum)return
//            Math.pow(1.0*n/maxMemoryNum/2,2.0)*4;//+0.5;
//            Math.pow(2,1.0*n/maxMemoryNum)/2;
            Math.max(0,1.0*n/maxMemoryNum-0.5)+0.5;
//            1.0*n/maxMemoryNum;
//            0.5*n/maxMemoryNum+0.5;
//            Math.log(n)/Math.log(maxMemoryNum);
//            1;
        int[] capacity = KLLSketchLazyExact.calcLevelMaxSize(maxMemoryNum,n);
        int level = capacity.length;
        int K = capacity[level-1];
        double CDF_COEF=2.296,CDF_EXP=0.9723;
        double esti_err = CDF_COEF / Math.pow(K, CDF_EXP);
//            System.out.println("n:"+n+"\t\t\t\t\tlv:"+level+"\t\testi_rank_err:"+(int)(n*esti_err));
        int maxERR = 0;
        int LazyNum = Math.max(0,capacity[level-1]-(n>>>(level-1)));
        for(int i=0;i<level-1;i++)maxERR+=(n/((capacity[i] + 1+LazyNum) <<i))<<i;

        int nextN=(int)(n*esti_err*2);
//            if(maxERR*2<nextN)System.out.println("\t\t\t\t\terr estimate...  n:"+n+"\t\tmemNum:"+maxMemoryNum+"\t\t\tbaselineEsti:"+nextN+"\t\t1.0Bound:"+maxERR*2+"\t\t[0,lv):"+level+"\t\tLazyNum:"+LazyNum);
        nextN=Math.min(nextN,maxERR*2);
        return 1.0+estimateIteration(nextN,maxMemoryNum);
    }

//    private double simulateIteration(double casePr,int depth,int n, int maxMemoryByte){
//        int maxMemoryNum = maxMemoryByte/8;
//        if(n<=0)return 0;
//        if(n<=maxMemoryNum)return
//            Math.max(0.8,Math.log(n)/Math.log(maxMemoryNum));
////            Math.max(0,1.0*n/maxMemoryNum-0.5)+0.5;
////            1;
//        int[] compactNum = new KLLSketchLazyEmptyForSimuCompact(n,maxMemoryNum).simulateCompactNumGivenN(n);
//        int maxERR=0;
//        for(int i=0;i<compactNum.length;i++)maxERR+=compactNum[i]<<i;
//        double bestSimuIter=1e3,bestSimuPr=1.0;
//
//        if(casePr>=0.01) {
//            ObjectArrayList<DoubleIntPair> prErrList=new ObjectArrayList<>();
//            int lastPrErr=Integer.MAX_VALUE;
//            for (double pr : simuPrList){
//                int prERR = KLLSketchLazyExact.queryRankErrBound(compactNum, pr);
//                if(prERR==lastPrErr/*||(prERR<=maxMemoryNum&&lastPrErr<=maxMemoryNum)*/)prErrList.set(prErrList.size()-1,DoubleIntPair.of(pr,prERR));
//                else {
//                    prErrList.add(DoubleIntPair.of(pr, prERR));
//                    lastPrErr = prERR;
//                }
//            }
//            for (DoubleIntPair prERRPair : prErrList) {
//                double Pr=prERRPair.leftDouble();int prERR=prERRPair.rightInt();
//                int nextN = prERR * 2, failN = (Math.min(n, maxERR * 2) - prERR) / 2;
//                double cntSimuIter = 1 + Pr * simulateIteration(casePr * Pr,depth+1, nextN, maxMemoryByte) + (1 - Pr) * (1 + simulateIteration(0*casePr * (1 - Pr), depth+1,failN, maxMemoryByte));
//                if (cntSimuIter <= bestSimuIter) {
//                    bestSimuIter = cntSimuIter;
//                    bestSimuPr = Pr;
//                }
//            }
////            System.out.println("\t\t\tsimulating  n:" + n + "\tcasePr:"+casePr+"\tdepth:"+depth+"\tbestSimuPR:" + bestSimuPr + "\tbestSimuIter:" + bestSimuIter+"\t\tsizeOfALLSimuPrList:"+simuPrList.size()+"\t\ttriedPrList:"+prErrList.size()+"\t"+prErrList);
//        }else{
//            double Pr = 0.99;
//            int prERR = KLLSketchLazyExact.queryRankErrBound(compactNum, Pr);
//            int nextN = prERR * 2, failN = (Math.min(n, maxERR * 2) - prERR) / 2;
//            bestSimuIter = 1 + Pr * simulateIteration(casePr * Pr,depth+1, nextN, maxMemoryByte) + (1 - Pr) * (1 + simulateIteration(0*casePr * (1 - Pr), depth+1,failN, maxMemoryByte));
////                    System.out.println("\t\t\tsimulating FAILED  n:"+n+"\tdepth:"+depth+"\tfixPR:"+Pr+"\tprERR:"+prERR+"\t\tsucc,failN:"+nextN+","+failN+"\t\tsimuIter:"+bestSimuIter);
//        }
//        return bestSimuIter;
//    }
    KLLSketchLazyEmptyForSimuCompact simuWorker;
    private double[] simulateIteration(double casePr,double lastPr,int depth,int n, int maxMemoryNum,int[] compactNum) {
        if (n <= 0) return new double[]{0,1.0};
        if (n <= maxMemoryNum) return
            new double[]{
            Math.max(0.75, Math.log(n) / Math.log(maxMemoryNum)),
                1.0};
        int maxERR = 0;
        for (int i = 0; i < compactNum.length; i++) maxERR += compactNum[i] << i;
        double bestSimuIter = 1e3, bestSimuPr = 1.0;

        if (casePr >= 0.01) {
            double[] prList = new double[simuPrList.size()];
            int[] prErrList = new int[simuPrList.size()];
            int lastPrErr = Integer.MAX_VALUE,prCount=0;
//            ObjectArrayList<DoubleIntPair> prErrList = new ObjectArrayList<>(simuPrList.size());
            for (double pr : simuPrList) if(pr>=lastPr||(lastPr==1.0&&pr==simuPrList.getDouble(simuPrList.size()-1))){
                int prERR = KLLSketchLazyExact.queryRankErrBound(compactNum, pr);
                if (prERR == lastPrErr)
                    prList[prCount-1]=pr;
//                    prErrList.set(prErrList.size() - 1, DoubleIntPair.of(pr, prERR));
                else {
//                    prErrList.add(DoubleIntPair.of(pr, prERR));
                    prList[prCount]=pr;
                    prErrList[prCount]=prERR;
                    prCount++;
                    lastPrErr = prERR;
                }
            }
            int[] successN = new int[prCount],failN = new int[prCount];
            int[][] succComNum= new int[prCount][],failComNum=new int[prCount][];
            for(int i=0;i<prCount;i++){
                int prERR = prErrList[i];
                successN[i]=prERR * 2+(prERR+prERR)/16;
                if(successN[i]<=maxMemoryNum)successN[i]+=(prERR+prERR)/16;
                failN[i] = (Math.min(n, maxERR * 2) - prERR) / 2;
            }
//            KLLSketchLazyEmptyForSimuCompact simuWorker = new KLLSketchLazyEmptyForSimuCompact(/*n, */maxMemoryNum);
            for(int i=0;i<prCount;i++)
                succComNum[i]=simuWorker.simulateCompactNumGivenN(successN[i]);
            for(int i=prCount-1;i>=0;i--)
                failComNum[i]=simuWorker.simulateCompactNumGivenN(failN[i]);
            for(int i=0;i<prCount;i++){
                double pr = prList[i];
                double cntSimuIter = 1 + pr * simulateIteration(casePr * pr, pr, depth + 1, successN[i], maxMemoryNum,succComNum[i])[0] + (1 - pr) * (1 + simulateIteration(0 * casePr * (1 - pr), pr,depth + 1, failN[i], maxMemoryNum,failComNum[i])[0]);
                if (cntSimuIter <= bestSimuIter) {
                    bestSimuIter = cntSimuIter;
                    bestSimuPr = pr;
                }
            }
//            if(depth==1)
//                System.out.println("\t\t\t\t\t\t\t\t\t\t\t\t//simulating  n:" + n + "\tcasePr:"+casePr+"\tdepth:"+depth+"\tbestSimuPR:" + bestSimuPr + "\tbestSimuIter:" + bestSimuIter+"\t\tsizeOfALLSimuPrList:"+simuPrList.size()+"\t\ttriedPrCount:"+prCount+"\t\tprERRUnderBestSimuPr:"+KLLSketchLazyExact.queryRankErrBound(new KLLSketchLazyEmptyForSimuCompact(maxMemoryNum).simulateCompactNumGivenN(n),bestSimuPr));
        }else{
            double pr = bestSimuPr = 0.99;
            int prERR = KLLSketchLazyExact.queryRankErrBound(compactNum, pr);
            int successN=prERR * 2;
            int failN = (Math.min(n, maxERR * 2) - prERR) / 2;
//            KLLSketchLazyEmptyForSimuCompact simuWorker = new KLLSketchLazyEmptyForSimuCompact(/*n, */maxMemoryNum);
            int[] succComNum,failComNum;
            if(failN<successN){
                failComNum = simuWorker.simulateCompactNumGivenN(failN);
                succComNum = simuWorker.simulateCompactNumGivenN(successN);
            }else{
                succComNum = simuWorker.simulateCompactNumGivenN(successN);
                failComNum = simuWorker.simulateCompactNumGivenN(failN);
            }
            bestSimuIter = 1 + pr * simulateIteration(casePr * pr, pr, depth + 1, successN, maxMemoryNum,succComNum)[0] + (1 - pr) * (1 + simulateIteration(0 * casePr * (1 - pr), pr, depth + 1, failN, maxMemoryNum,failComNum)[0]);
        }
        DEBUG_COUNT_SIMULATION++;
        return new double[]{bestSimuIter,bestSimuPr};
    }


//    private double evaluatePr(KLLSketchLazyExact sketch,int maxMemoryNum,long rk1,long rk2,double Pr,double[] deterministic_result){
//        double[] cntResult = sketch.findResultRange(rk1,rk2,Pr);
//        int rkValL = (int)cntResult[2],rkValR=(int)cntResult[3],prErrL=(int)cntResult[4],prErrR=(int)cntResult[5];
//
////        int nextN = Math.max(rkValR-rkValL,prErrL);
////        int nextN = (prErrL+prErrR);
//        int nextN = Math.max(rkValR-rkValL,prErrL+prErrR);
//        nextN+=(prErrL+prErrR)/16;
//        if(nextN<=maxMemoryNum)
//            nextN+=(prErrL+prErrR)/32;
////        int nextN = (rkValR-rkValL+(prErrL+prErrR)/4);
//        int failN = (((int)deterministic_result[3]-(int)deterministic_result[2])-nextN)/2;
////        double estimateIterationNum = Pr*estimateIteration(nextN,maxMemoryNum)+(1-Pr)*(1+estimateIteration(failN,maxMemoryNum));
//        double simulateIterationNum = Pr*simulateIteration(Pr,1,nextN,maxMemoryNum)+(1-Pr)*(1+simulateIteration(0*(1-Pr),1,failN,maxMemoryNum));
//        System.out.println("\t\t\t\t\t\t\t\tcntPR:"+Pr+"\tsuccessN:\t"+nextN+"\tprERR:"+prErrL+"\t\tfailN:\t"+failN+/*"\t\testi_iter:\t"+estimateIterationNum+*/"\t\tsimu_iter:\t"+simulateIterationNum);
////        System.out.println("\t\t\t\testi_iter_num:\t"+estimateIterationNum);
////        return estimateIterationNum;
//        return simulateIterationNum;
//    }
    private double[] evaluatePr(int maxMemoryNum,double Pr,int succN,int failN,int[] succComNum,int[]failComNum){
        double[] simulateResult= new double[3];
        double[] successResult = simulateIteration(Pr,Pr,1,succN,maxMemoryNum,succComNum);
        double[] failResult = simulateIteration(0*(1-Pr),Pr,1,failN,maxMemoryNum,failComNum);
        simulateResult[0] = Pr*successResult[0]+(1-Pr)*(1+failResult[0]);
        simulateResult[1] = successResult[1];
        simulateResult[2] = failResult[1];
//        if(DEBUG_PRINT)System.out.println("\t\t\t\t\t\t\t\tcntPR:"+Pr+"\tsuccessN:\t"+succN+"\t\tfailN:\t"+failN+/*"\t\testi_iter:\t"+estimateIterationNum+*/"\t\tsimu_iter:\t"+simulateResult[0]+"\tsimu_nextSuccessPr:"+simulateResult[1]);
        return simulateResult;
    }
    private double findBestPr(KLLSketchLazyExact sketch, int maxMemoryNum,long rk1,long rk2,double[] deterministic_result){
        double bestPr=1.0,bestEstiNum=1e9,nextSuccessPr=1.0;
        int[] successN = new int[prList.size()],failN = new int[prList.size()];
        int[][] succComNum= new int[prList.size()][],failComNum=new int[prList.size()][];
        for(int i=0;i<prList.size();i++){
            double pr=prList.getDouble(i);
            double[] cntResult = sketch.findResultRange(rk1,rk2,pr);
            int rkValL = (int)cntResult[2],rkValR=(int)cntResult[3],prErrL=(int)cntResult[4],prErrR=(int)cntResult[5];
            int tmpSuccessN = Math.max(rkValR-rkValL,prErrL+prErrR);
            tmpSuccessN+=(prErrL+prErrR)/16;
            if(tmpSuccessN<=maxMemoryNum)
                tmpSuccessN+=(prErrL+prErrR)/16;
            int tmpFailN = (((int)deterministic_result[3]-(int)deterministic_result[2])-tmpSuccessN)/2;
            successN[i]=tmpSuccessN;
            failN[i]=Math.max(tmpSuccessN,tmpFailN);
        }
        //KLLSketchLazyEmptyForSimuCompact
        simuWorker = new KLLSketchLazyEmptyForSimuCompact(/*(int)sketch.getN()/2, */maxMemoryNum);
        for(int i=0;i<prList.size();i++)
            succComNum[i]=simuWorker.simulateCompactNumGivenN(successN[i]);
        for(int i=prList.size()-1;i>=0;i--)
            failComNum[i]=simuWorker.simulateCompactNumGivenN(failN[i]);

        for(int i=0;i<prList.size();i++) {
            double[] cntPrResult = evaluatePr(maxMemoryNum, prList.getDouble(i), successN[i],failN[i],succComNum[i],failComNum[i]);
            if(cntPrResult[0]<=bestEstiNum){
                bestEstiNum = cntPrResult[0];
                bestPr = prList.getDouble(i);
                nextSuccessPr = cntPrResult[1];
            }
//            System.out.println("\t\t\t\t\tcntPR:"+prList.getDouble(i)+"\tsuccessN:\t"+successN[i]+"\t\tfailN:\t"+failN[i]+"\t\tsimu_iter:\t"+cntNum);
        }
        if(DEBUG_PRINT)System.out.println("bestPr:"+bestPr+"\t\testiIter:"+bestEstiNum+"\t\t\tnextSuccessPr:\t"+nextSuccessPr);
        return bestPr;
    }


    public void testError(int queryN, int maxMemoryByte) throws IOException {
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

//        System.out.println("\t!! DataSketchesKLL_K="+DataSketchesKLL_K);
        for (int T = 0; T < TEST_CASE; T++) {
//            System.out.println("\tTEST_ID:"+T);
//            int L = 0, R = queryN;
            int L=LL[T],R=RR[T];

//            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
//            Arrays.sort(query_a);

            double q_tmp = 1e-2, q_start = q_tmp, q_end = 1 - q_tmp;
                q_start=0.50;q_end=0.50;
            double q_add = 1e-2, q_count = Math.floor((q_end - q_start + 1e-10) / q_add) + 1;
            double ratioPerQuery = 1.0 / (q_count * TEST_CASE);
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
//                if(q<0.189||q>0.194)continue;
//                int query_rank1 = (int) Math.floor(q * queryN),query_rank2 = (int) Math.ceil(q * queryN);
//                int query_rank1 = (int) Math.round(q * queryN),query_rank2 = query_rank1;
                int query_rank1 = (int) Math.floor(q * (queryN-1)+1),query_rank2 = (int) Math.ceil(q * (queryN-1)+1);
//                System.out.println("\n\t\t\t\t\t\tq rank1,2:"+q+"\t"+query_rank1+","+query_rank2);
                long last_n=queryN;


//                double[] deterministic_result,iterate_result = lazyWorker.findResultRange(query_rank1,query_rank2,tmpPR);
//                deterministic_result = lazyWorker.findResultRange(query_rank1,query_rank2,1.0);
//                avg_iteration+=ratioPerQuery;
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
                    for(int i=L;i<R;i++) {
                        if (a[i] >= valL && a[i] <= valR) cntWorker.update(dataToLong(a[i]));
                        else if (a[i] < valL) CountOfLessThanValL++;
                    }
//                    if(MMP==1&&T==0&&q==q_start)cntWorker.showCompact();
//                    if(cntWorker.getN()==cntWorker.getNumLen())System.out.println("\t\t\tfinal iter\tspaceRate:"+1.0*cntWorker.getN()/(maxMemoryByte/8));
//                    cntWorker.show();
//                    System.out.println("\t\t\t\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntSketch_N:"+cntWorker.getN());
                    int cntRank1=query_rank1-CountOfLessThanValL;
                    int cntRank2=query_rank2-CountOfLessThanValL;
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

                    double bestPr;
                    if(prList==null) {
//                        prList = new DoubleArrayList();
//                        for (double tmp = 0.5; tmp >/*1e-15*/1e-3 + 1e-6; tmp *= 0.85/*tmp/=2.0*/) {
//                            prList.add(1 - tmp);
//                        }
//                        prList.add(0.999);
                        prList = new DoubleArrayList();
                        for(double tmp=0.70;tmp<0.9-1e-6;tmp+=0.02)prList.add(tmp);
                        for(double tmp=0.90;tmp<0.99-1e-6;tmp+=0.01)prList.add(tmp);
                        prList.add(0.99);
                        prList.add(0.995);
                        prList.add(1.0);
//                        for(double tmp=0.99;tmp<0.999+1e-6;tmp+=1e-3)prList.add(tmp);
                    }
                    if(simuPrList==null){
                        simuPrList = new DoubleArrayList();
                        for(double pr:prList)
                            if(pr<1.0)
                                simuPrList.add(pr);
//                        for (double tmp = 0.5; tmp >/*1e-15*/1e-2 + 1e-6; tmp *= 0.85/*tmp/=2.0*/) {
//                            simuPrList.add(1 - tmp);
//                        }
//                        simuPrList.add(0.999);
                    }
                    bestPr = findBestPr(cntWorker,maxMemoryByte/8,cntRank1,cntRank2,deterministic_result);
//                    for(double tmp:prList) {
//                        double cntNum = evaluatePr(cntWorker, maxMemoryByte/8, cntRank1, cntRank2, tmp, deterministic_result);
//                        if(cntNum<=bestEstiNum){
//                            bestEstiNum = cntNum;
//                            bestPr = tmp;
//                        }
//                    }
////                    bestPr=0.95; // DANGER
                    if(DEBUG_PRINT)System.out.println("bestPr:"+bestPr+"\t\tcntIter:"+MMP);

                    ALLPrSum+=bestPr;ALLPrCount+=1;
                    iterate_result = cntWorker.findResultRange(cntRank1,cntRank2,bestPr);
//                    System.out.println("\t\t\t\tbestPr:"+bestPr+"\t\t\tcntN:"+cntWorker.getN());
//                    System.out.println("\t\t\t\t\tcntL,R:"+iterate_result[0]+","+iterate_result[1]+"\t\t\tCountOfLessThanValL:"+CountOfLessThanValL+"\t\tcntRank:"+cntRank1+" "+cntRank2);
                }
                double exact_quantile_v = (deterministic_result[0]+deterministic_result[1])*0.5;
                if(DEBUG_PRINT)System.out.println("FINISH CALC\t\tq:"+q+"\texact_quantile="+exact_quantile_v+"\t\titer:\t"+MMP+"\n");
//                if(Math.abs(exact_quantile_v-N*q)>1.0){
//                    System.out.println("!!!!!!!!!!ERROR!!!!");
//                    return;
//                }
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("DynamicPR!\tTEST_CASE="+TEST_CASE+"\tDATASET:"+dataType+"\tqueryN:\t" + queryN+"\tmemory:\t" + maxMemoryByte+"\t\tavg_iteration:\t"+avg_iteration+"\t\tavgPrChose:"+ALLPrSum/ALLPrCount+"\t\tavgFailRate:"+(FailCount/IterCount)+"\tusedNorDistri:\t"+KLLSketchLazyExact.DEBUG_COUNT_USENORMAL+"\tnot_cached_calc_pr_err:"+KLLSketchLazyExact.DEBUG_COUNT_USEPRERR+"\t\tSIMU_COUNT:"+DEBUG_COUNT_SIMULATION);
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START_T = new Date().getTime();
        MainForExactLazyKLLFindPr main;

        main = new MainForExactLazyKLLFindPr();
        for (int dataType = startType; dataType <= endType; dataType++) { // CHECK IT
            main.RESULT_LINE=0;
            main.prepareA(dataType);
            for (int queryN : new int[]{/*10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000,50000000/*,1000000,1500000,2000000,2500000,*/3000000,3500000,4000000,4500000,5000000/*,5500000,6000000/**/})
//            for (int queryN : new int[]{3000000/*,100000,200000,400000,500000,800000,1000000,2000000,3000000/*,5000000,10000000,20000000,50000000/**/})
                for (int query_mem : new int[]{/*1024*4*//*1024*1,1024*2,*/1024*4/*,1024*8,1024*16,1024*32/**/})
                    main.testError(queryN, query_mem);
//                    for(double tmpPr=0.001;tmpPr>=1e-10;tmpPr/=2) {
//                        main.tmpPR=1-tmpPr;
//                        main.testError(queryN, query_mem);
//                    }
                main.RESULT_LINE++;
//            System.out.println("byPage & byChunkDivide\nTEST_CASE=" + TEST_CASE);
        }
//        System.out.println("\nError rate:");
//        for (String s : main.err_result)
//            System.out.println(s);
        System.out.println("\t\t\tALL_TIME:"+(new Date().getTime()-START_T));
    }
}

import com.carrotsearch.hppc.ObjectStack;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.eclipse.collections.api.tuple.primitive.IntLongPair;
import org.eclipse.collections.impl.tuple.primitive.PrimitiveTuples;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class IntervalEvaluatingKLL_LSM {
    int dataType=-233;
    static int baseSketchCompression = 32,LSM_WIDTH=30,addKLLLevel=3;
    static int pagePerChunkL=1,pagePerChunkR=1;
    static int chunkN/* = 4096*16 */, chunkNum/* = N/chunkN*/, pageN = 4096;
    static int N = 81000000/pageN*pageN, pageNum = N / pageN; // CHECK IT  Can't query 5E7.
    public static int TEST_CASE = 512; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] KLLWorkerByChunkDivide;
    static ObjectArrayList<ObjectArrayList<KLLSketchForQuantile>> KLLWorkerForLSMNode;
    ArrayList<String> err_result = new ArrayList<>();
    ArrayList<String> time_result = new ArrayList<>();
    boolean show_time = false, show_err = true;
    int RESULT_LINE = 0;
    Random random = new Random(233);

    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = //longToResult(i&4095);
                        //Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
                    random.nextDouble();
            return;
        }
        BufferedReader reader = null;
        if (dataType == 1)
            reader = new BufferedReader(new FileReader(new File("1_bitcoin.csv")));
        if (dataType == 2)
            reader = new BufferedReader(new FileReader(new File("SpacecraftThruster.txt")));
        if (dataType == 3)
            reader = new BufferedReader(new FileReader(new File("4_taxipredition8M.txt")));
        if (dataType == 4)
            reader = new BufferedReader(new FileReader(new File("5_wh.csv")));
        assert reader != null;
        reader.readLine(); // ignore first line.
        String line;
        int cntN = 0;
        while ((line = reader.readLine()) != null) {
            a[cntN++] = Double.parseDouble(line);
            if (cntN == N) break;
        }
    }


    public void prepareWorker(int maxSeriesByte) {
        int pageNumInChunk =chunkN/pageN;
        int compactionLevel = 0;
        for(int i=0;i<=10;i++)if((maxSeriesByte<<i)>=pageN*8){
            compactionLevel=i;
            break;
        }
//        KLLWorkerByPage = new KLLSketchForQuantile[pageNum];
        KLLWorkerByChunkDivide= new KLLSketchForQuantile[pageNum];
        KLLWorkerForLSMNode = new ObjectArrayList<>();
        for(int i=0;i<chunkNum;i++){
            ObjectArrayList<IntLongPair> TVList = new ObjectArrayList<>();
            for(int j=i*chunkN;j<i*chunkN+chunkN;j++){
                int index = j;
                TVList.add(PrimitiveTuples.pair(index,dataToLong(a[index])));
            }
            TVList.sort(Comparator.comparingInt(IntLongPair::getOne));
            LongArrayList t=new LongArrayList(chunkN),v=new LongArrayList((chunkN));
            for(int j=0;j<chunkN;j++){
                t.add(TVList.get(j).getOne());
                v.add(TVList.get(j).getTwo());
            }
            IntArrayList pageStartIndex = new IntArrayList(pageNumInChunk);
            LongArrayList pageMinV = new LongArrayList(pageNumInChunk),pageMaxV = new LongArrayList(pageNumInChunk);
            for(int j=0;j<pageNumInChunk;j++) {
                pageStartIndex.add(j * pageN);
                long minV=Long.MAX_VALUE,maxV=Long.MIN_VALUE;
                for(int k=j*pageN;k<j*pageN+pageN;k++){
                    minV=Math.min(minV,v.getLong(k));
                    maxV=Math.max(maxV,v.getLong(k));
                }
                pageMinV.add(minV);
                pageMaxV.add(maxV);
            }
            KLLSketchForDivide chunkSketch = new KLLSketchForDivide(compactionLevel,t,v);
//            for(int j=i*chunkN/pageN;j<i*chunkN/pageN+chunkN/pageN;j++)
//                KLLWorkerByChunkDivide[j] = new LongKLLSketchForUnit(chunkSketch,workerMaxT[j],pageN);
            List<LongKLLSketchForUnit> pageSketch = chunkSketch.divideMemSketchByItemValue(pageStartIndex,pageMinV,pageMaxV);
            for(int j=0;j<pageNumInChunk;j++) {
                KLLWorkerByChunkDivide[i * pageNumInChunk + j] = pageSketch.get(j);
//                pageSketch.get(j).showNum();
            }
//            if(i==0)KLLWorkerByChunkDivide[0].showNum();
        }
//        for(int i=0;i<100;i++)System.out.println(workerMinT[i]+"---"+workerMaxT[i]);

        KLLWorkerForLSMNode.add(new ObjectArrayList<>());


        for(int i=0;i<pageNum;i++){
            KLLWorkerForLSMNode.get(0).add(KLLWorkerByChunkDivide[i]);
            int lvSize = KLLWorkerForLSMNode.get(0).size();
            if(lvSize%(LSM_WIDTH*pageNumInChunk)==0)
                for(int lv=0;lv<KLLWorkerForLSMNode.size()&&lvSize%(LSM_WIDTH)==0;lv++){
                    if(lv+1==KLLWorkerForLSMNode.size())
                        KLLWorkerForLSMNode.add(new ObjectArrayList<>());
                    KLLSketchForLSMFile nextSketch = new KLLSketchForLSMFile();
                    int sketchNum = lv==0?(LSM_WIDTH*pageNumInChunk):LSM_WIDTH;
                    for(int k=0;k<sketchNum;k++)
                        nextSketch.addSubSketch(KLLWorkerForLSMNode.get(lv).get(lvSize-sketchNum+k));
                    nextSketch.compactSubSketches(addKLLLevel);
                    KLLWorkerForLSMNode.get(lv+1).add(nextSketch);
                    lvSize = KLLWorkerForLSMNode.get(lv+1).size();
                }
        }
//        for(int lv=1;lv<KLLWorkerForLSMNode.size();lv++){
//            System.out.println("\t\tLSM level "+lv+"\t\tlvNodeNum:"+KLLWorkerForLSMNode.get(lv).size());
//            long lvN=0;
//            ObjectArrayList<KLLSketchForQuantile> nodes = KLLWorkerForLSMNode.get(lv);
//            for(KLLSketchForQuantile node :nodes)
//                lvN+=node.getN();
//            System.out.println("\t\t\tlvN:"+lvN);
//            for(KLLSketchForQuantile node :nodes) {
//                node.show();break;
////                node.showNum();
//            }
//            System.out.println("------------------------------------");
//        }
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

    private boolean overlapInterval(long x, long y, long L, long R) { // [L,R]
        return !inInterval(x, y, L, R) && !(y < L || x > R);
    }


    public void testError(int queryN, int maxMemoryByte, int maxSeriesByte) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        double err_full = 0, err_merge_lsm = 0, err_merge_divide = 0;
        double computation_full = 0, computation_merge_lsm = 0, computation_merge_divide = 0;
        double[] query_a = new double[queryN];
        int[] LL = new int[TEST_CASE];
        int[] RR = new int[TEST_CASE];
        Random random = new Random(233);
        for (int i = 0; i < TEST_CASE; i++) {
            LL[i] = random.nextInt(N - queryN + 1);
            RR[i] = LL[i] + queryN;
        }
        DoubleArrayList errLSM=new DoubleArrayList(),errDivide=new DoubleArrayList(),errFull=new DoubleArrayList();

        int pagePerChunk = chunkN/pageN;
        for (int T = 0; T < TEST_CASE; T++) {
            if((T&1)==0)//CHECK IT!
                prepareWorker(maxSeriesByte);
            int L = LL[T], R = RR[T];
            int pageL = (L+pageN-1)/pageN, pageR = R/pageN;
            int pagePosL = pageL*pageN, pagePosR = pageR*pageN;
//            int chunkL = (L+chunkN-1)/chunkN, chunkR = R/chunkN-1;
//            int chunkPosL = chunkL*chunkN, chunkPosR = chunkR*chunkN;
//            System.out.println("\t\t\t"+pagePosL+"\t"+pagePosR);

            HeapLongStrictKLLSketch merge_lsm_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for(int i=L;i<Math.min(R,pagePosL);i++) {
                merge_lsm_worker.update(dataToLong(a[i]));
                computation_merge_lsm++;
            }
            ObjectArrayList<KLLSketchForQuantile> tmpList1=new ObjectArrayList<>();
            ObjectStack<KLLSketchForQuantile>tmpList2=new ObjectStack<>();
            for(int tmpL=pageL,tmpR=pageR-1,lv=0;tmpL<=tmpR;lv++){
                int lvWidth = lv==0?LSM_WIDTH*pagePerChunk:LSM_WIDTH;
//                System.out.println("\t\t\ttmpLR:"+tmpL+"\t"+tmpR);
                if(tmpL/lvWidth==tmpR/lvWidth) {
                    if(tmpR-tmpL+1==LSM_WIDTH)
                        tmpList1.add(KLLWorkerForLSMNode.get(lv+1).get(tmpL/LSM_WIDTH));
                    else
                    for(int i=tmpL;i<=tmpR;i++)
                        tmpList1.add(KLLWorkerForLSMNode.get(lv).get(i));
                    break;
                }
                if(tmpL%lvWidth>0)
                for(int i=0;i<lvWidth-tmpL%lvWidth;i++)
                    tmpList1.add(KLLWorkerForLSMNode.get(lv).get(tmpL+i));
//                System.out.println("\t\t\t\tlist1Size="+tmpList1.size());

                if(tmpR%lvWidth<lvWidth-1)
                for(int i=0;i<=tmpR%lvWidth;i++) {
                    tmpList2.push(KLLWorkerForLSMNode.get(lv).get(tmpR-i));
                }
//                System.out.println("\t\t\t\tlist2Size="+tmpList2.size());
                tmpL=(tmpL+lvWidth-1)/lvWidth;
                tmpR=(tmpR-lvWidth+1)/lvWidth;
            }
//            System.out.println("\t\t\t\tlist1Size="+tmpList1.size());
//            System.out.println("\t\t\t\tlist2Size="+tmpList2.size());
            for(KLLSketchForQuantile sketch:tmpList1) {
                merge_lsm_worker.mergeWithTempSpace(sketch);
                computation_merge_lsm+=sketch.getNumLen();
//                System.out.println("\t\t\t\t\t\t"+sketch.getNumLen());
            }
            while(!tmpList2.isEmpty()) {
                KLLSketchForQuantile sketch = tmpList2.pop();
                merge_lsm_worker.mergeWithTempSpace(sketch);
                computation_merge_lsm+=sketch.getNumLen();
//                System.out.println("\t\t\t\t\t\t"+sketch.getNumLen());
            }
            for(int i=Math.max(L,pagePosR);i<R;i++) {
                merge_lsm_worker.update(dataToLong(a[i]));
                computation_merge_lsm++;
            }
//            System.out.println("\t\tquery L,R:"+L+","+R+"\t\tmerged_lsm_N:"+merge_lsm_worker.getN());

            HeapLongStrictKLLSketch merge_divide_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for(int i=L;i<Math.min(R,pagePosL);computation_merge_divide++,i++)
                merge_divide_worker.update(dataToLong(a[i]));
            for (int i=pageL;i<pageR; i++) {
                merge_divide_worker.mergeWithTempSpace(Collections.singletonList(KLLWorkerByChunkDivide[i]));
                computation_merge_divide+=KLLWorkerByChunkDivide[i].getNumLen();
            }
            for(int i=Math.max(L,pagePosR);i<R;computation_merge_divide++,i++)
                merge_divide_worker.update(dataToLong(a[i]));

            HeapLongStrictKLLSketch full_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
            for(int i=L;i<R;computation_full++,i++)
                full_worker.update(dataToLong(a[i]));
//            System.out.println("\t\tQuerySketchNumLen:"+merge_lsm_worker.getNumLen()+"\t"+merge_divide_worker.getNumLen()+"\t"+full_worker.getNumLen());
            if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
            Arrays.sort(query_a);


//            merge_lsm_worker.show();
//            merge_divide_worker.show();

            double q_start = 0.01, q_end = 0.99, q_add = 0.0001, q_count = Math.floor((q_end - q_start - 1e-10) / q_add) + 1;
            double ratioPerQuery = 1.0 / (q_count * TEST_CASE);
            for (double q = q_start; q < q_end + 1e-10; q += q_add) {
                int query_rank = (int) (q * queryN);

                double merge_lsm_v = longToResult(merge_lsm_worker.findMinValueWithRank(query_rank));
                int merge_lsm_delta_rank = getValueActualRank(query_a, queryN, merge_lsm_v) - query_rank;
                double merge_lsm_relative_err = 1.0 * merge_lsm_delta_rank / (queryN);
                err_merge_lsm += Math.abs(merge_lsm_relative_err) *ratioPerQuery;

                double merge_divide_v = longToResult(merge_divide_worker.findMinValueWithRank(query_rank));
                int merge_divide_delta_rank = getValueActualRank(query_a, queryN, merge_divide_v) - query_rank;
                double merge_divide_relative_err = 1.0 * merge_divide_delta_rank / (queryN);
                err_merge_divide += Math.abs(merge_divide_relative_err) *ratioPerQuery;

                double full_v = longToResult(full_worker.findMinValueWithRank(query_rank));
                int full_delta_rank = getValueActualRank(query_a, queryN, full_v) - query_rank;
                double full_relative_err = 1.0 * full_delta_rank / (queryN);
                err_full += Math.abs(full_relative_err) *ratioPerQuery;

                errLSM.add(Math.abs(merge_lsm_relative_err));
                errDivide.add(Math.abs(merge_divide_relative_err));
                errFull.add(Math.abs(full_relative_err));
//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
//                    System.out.println("\t\t\trk:"+query_rank+"\t\t"+Math.abs(merge_divide_relative_err));
            }

        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("\t\t" + queryN+","+maxMemoryByte + "\t" + err_merge_lsm + "\t" + err_merge_divide + "\t" + err_full);
        //        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));
        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t\t\t\t\t" + err_merge_lsm + "\t" + err_merge_divide + "\t" + err_full));

        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t|\t" + calcSTDDEV(errLSM) + "\t" + calcSTDDEV(errDivide) + "\t" + calcSTDDEV(errFull)+"\t|\t"));
        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t|\t" + computation_merge_lsm/TEST_CASE + "\t" + computation_merge_divide/TEST_CASE + "\t" + computation_full/TEST_CASE+"\t|\t"));
    }

    private double calcSTDDEV(DoubleArrayList err){
        double avg = 0;
        for(double x:err)avg+=x;
        avg/=err.size();
        double STDDEV=0;

        for(double x:err)STDDEV+=(x-avg)*(x-avg);
        STDDEV/=err.size();
        STDDEV = Math.sqrt(STDDEV);
        return STDDEV;
    }


    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        long START=new Date().getTime();
        IntervalEvaluatingKLL_LSM main;

        System.out.println("range query: lsm & dividedPageSketch & fullCalc" + "\n");
        System.out.println("\tLSM_WIDTH="+LSM_WIDTH+"\taddKLLLevel="+addKLLLevel);
        for (int dataType = 2; dataType <= 2; dataType++) { // CHECK IT
            main = new IntervalEvaluatingKLL_LSM();
            main.prepareA(dataType);

            assert pagePerChunkL==pagePerChunkR;
            int PagePerChunk=pagePerChunkL;
//            for (int queryN : new int[]{1000000,2000000,4000000,10000000,20000000,40000000})
            for (int queryN : new int[]{5000000,10000000,20000000,40000000,80000000})
                for (int query_mem : new int[]{/*1024*32,1024*64,*/1024*128/*,1024*256,1024*512*/}){
                chunkN=pageN*PagePerChunk;
                chunkNum = N/chunkN;
                main.err_result.add("query_mem:\t" + query_mem + "\tqueryN:\t"+queryN);
                main.testError(queryN, query_mem, pageN*8/baseSketchCompression);
                main.RESULT_LINE++;
            }
            System.out.println("lsm & dividedPageSketch & fullCalc\nTEST_CASE=" + TEST_CASE+"\tpageN="+pageN+"\tPagePerChunk="+PagePerChunk+"\tsketchCompression="+baseSketchCompression+"\tdataset:"+dataType+"\tLSM_WIDTH="+LSM_WIDTH);
            System.out.println("\nError rate:           err_lsm    err_divideChunk   err_NoPreComputation |   theirSTDDEV");
            for (String s : main.err_result)
                System.out.println(s);
            System.out.println("\n------------------------------------------------------------------\n\n");
        }
        System.out.println("\t\tALL_TIME:"+(new Date().getTime()-START));
    }
}

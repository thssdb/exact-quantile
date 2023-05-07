import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

/** Compare sketchByPage & sketchByChunkDivide.
 * Different ChunkSize.
*/
public class MainForKLLSST {
    int dataType;
    static int LSM_T=30;
    static int pageN = 8192;
    static int N = 55000000/pageN*pageN, pageNum=N/pageN; // CHECK IT
    public static int TEST_CASE = 1; // CHECK IT
    static double[] a;
    static KLLSketchForQuantile[] pageKLL;
    static ObjectArrayList<ObjectArrayList<SST>> lsmNode;
    static long[] workerMinT, workerMaxT;
    ArrayList<String> err_result = new ArrayList<>();
    ArrayList<String> time_result = new ArrayList<>();
    boolean show_time = false, show_err = true;
    int RESULT_LINE = 0;
    Random random = new Random(233);

    public class SST{
        public int level;
        public ObjectArrayList<ObjectArrayList<KLLSketchForQuantile>> sketch;
        public ObjectArrayList<LongArrayList> nodeMinT,nodeMaxT;
        public ObjectArrayList<IntArrayList>nodeStartDataID;
    }


    public void prepareA(int dataType) throws IOException {
        if (a == null) a = new double[N];
        this.dataType = dataType;

        if (dataType == 0) {
            for (int i = 0; i < N; i++)
                a[i] = longToResult(i);
                        //Math.pow(-1, random.nextInt(2)) * Math.pow(10.0, (2 * Math.pow(random.nextDouble(), 2) - 1) * 300);
//                    random.nextDouble();
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

    public void prepareWorker(int maxSeriesByte,int sketchSizeRatio){
        pageKLL = new KLLSketchForQuantile[pageNum];
        int enoughMemByte = pageN*10;
        for(int i=0;i<pageNum;i++) {
            LongKLLSketch worker = new LongKLLSketch(pageN, enoughMemByte, maxSeriesByte);
            for (int j = 0; j<pageN; j++) worker.update(dataToLong(a[i*pageN+j]));
            worker.compactBeforeSerialization();
            pageKLL[i] = worker;
//            if(i==0)worker.show();
        }

        int sketchNum = 1,SSTlv=0;
        lsmNode = new ObjectArrayList<>();
        lsmNode.add(new ObjectArrayList<>());
        while(sketchNum*LSM_T<pageNum) {
            sketchNum*=LSM_T;
            SSTlv++;
            lsmNode.add(new ObjectArrayList<>());
        }
//        System.out.println("\t\t??\t"+lsmNode.size());
        for(int remainingPage=pageNum;SSTlv>0;remainingPage%=LSM_T,SSTlv--,sketchNum/=LSM_T){
            for(int i=0;i<remainingPage/sketchNum;i++){
                SST sst = new SST();
                lsmNode.get(SSTlv).add(sst);
                sst.sketch = new ObjectArrayList<>();
                sst.nodeMinT = new ObjectArrayList<>();
                sst.nodeMaxT = new ObjectArrayList<>();
                sst.nodeStartDataID = new ObjectArrayList<>();

                for(int sketchLV=0,width=sketchNum;sketchLV<=SSTlv;sketchLV++,width/=LSM_T) {
                    sst.sketch.add(new ObjectArrayList<>());
                    sst.nodeMinT.add(new LongArrayList());
                    sst.nodeMaxT.add(new LongArrayList());
                    sst.nodeStartDataID.add(new IntArrayList());
                    if(sketchLV==0) {
                        for (int j = 0; j < width; j++) {
                            int pageID = pageNum - remainingPage + i * sketchNum + j;
                            sst.sketch.get(0).add(pageKLL[pageID]);
                            sst.nodeMinT.get(0).add(pageID * pageN);
                            sst.nodeMaxT.get(0).add(pageID * pageN + pageN - 1);
                            sst.nodeStartDataID.get(0).add(pageID * pageN);
                        }
                    }else{
                        for(int j=0;j<width;j++){
                            KLLSketchForSST sketch = new KLLSketchForSST();
                            for(int k=0;k<LSM_T;k++) {
                                sketch.addSubSketch(sst.sketch.get(sketchLV - 1).get(j * LSM_T + k));
                            }
                            sketch.compactSubSketches(sketchSizeRatio);
                            sst.sketch.get(sketchLV).add(sketch);
                            long mn = sst.nodeMinT.get(sketchLV-1).getLong(j * LSM_T);
                            long mx = sst.nodeMaxT.get(sketchLV-1).getLong(j * LSM_T+LSM_T-1);
                            sst.nodeMinT.get(sketchLV).add(mn);
                            sst.nodeMaxT.get(sketchLV).add(mx);
                            sst.nodeStartDataID.get(sketchLV).add(sst.nodeStartDataID.get(sketchLV-1).getInt(j * LSM_T));
//                            System.out.println("\t\ta sketch in sst.  T:"+mn+"..."+mx+"\t\t");
//                            System.out.print("\t\t\t\t");sketch.show();
                        }
                    }
                }
            }
        }
        for(int lv=lsmNode.size()-1;lv>=0;lv--)lsmNode.get(lsmNode.size()-1).get(0).sketch.get(lv).get(0).show();
    }


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

    private boolean inInterval(long x, long y, long L, long R) {
        return x >= L && y <= R;
    }

    private boolean inInterval(long x, long L, long R) {
        return x >= L && x <= R;
    }

    private boolean overlapInterval(long x, long y, long L, long R) { // [L,R]
        return !(y < L || x > R);
    }

    private void range_query(int lv,int p,SST sst,HeapLongStrictKLLSketch query_sketch,long L,long R,long otherL,long otherR){ // []

        long cntL=sst.nodeMinT.get(lv).getLong(p),cntR=sst.nodeMaxT.get(lv).getLong(p);
        if(!overlapInterval(cntL,cntR,L,R))return;
//            System.out.println("\t\tcntSST T:"+cntL+"..."+cntR);
        if(inInterval(cntL,cntR,L,R)&&inInterval(cntL,cntR,otherL,otherR)){
//            System.out.println("\t\tmerge with T:"+cntL+"..."+cntR);
            query_sketch.mergeWithTempSpace(Collections.singletonList(sst.sketch.get(lv).get(p)));
            return;
        }
        if(lv>0) {
            for(int i=0;i<LSM_T;i++)
            range_query(lv - 1, p*LSM_T+i, sst, query_sketch, L, R, otherL, otherR);
        }else{
            for (int start=sst.nodeStartDataID.get(lv).getInt(p),j = start; j < start + pageN; j++) {
                int index = j;//cc.getInt(j);
                if (inInterval(index, L, R))
                    query_sketch.update(dataToLong(a[index]));
            }
        }
    }

    public void testError(int queryN, int maxMemoryByte, int maxSeriesByte, int sketchSizeRatio) throws IOException {
//        System.out.println("\tdata/mem:"+queryPageNum*pageN*8/maxMemoryByte+
//            "\tpageData/pageKLL:"+pageN*8/maxSeriesByte);
        DecimalFormat fnum = new DecimalFormat("#0.00");
        double err_full = 0, err_merge_sst=0;
        double[] query_a = new double[queryN];


        for (int T = 0; T < TEST_CASE; T++) {
            prepareWorker(maxSeriesByte,sketchSizeRatio);

            for (int L = 0; L + queryN <= N; L += queryN) {
                int R = L+queryN;
//            System.out.println("\t\t\t"+posL+"\t"+posR);


                HeapLongStrictKLLSketch merge_sst_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
                for(int lsmLV=lsmNode.size()-1;lsmLV>=0;lsmLV--){
                    for(SST sst:lsmNode.get(lsmLV))
                        range_query(lsmLV,0,sst,merge_sst_worker,L,R-1,Long.MIN_VALUE,Long.MAX_VALUE);
                }

                HeapLongStrictKLLSketch full_worker = new HeapLongStrictKLLSketch(maxMemoryByte);
                for (int i = L; i < R; i++)
                    full_worker.update(dataToLong(a[i]));

                assert merge_sst_worker.getN()==full_worker.getN();

                if (R - L >= 0) System.arraycopy(a, L, query_a, 0, R - L);
                Arrays.sort(query_a);

                double q_add=0.0001,q_start=q_add,q_end=1-q_add,q_count = Math.floor((q_end-q_start-1e-10)/q_add)+1;
                double ratioPerQuery = 1.0/(q_count * TEST_CASE*N/queryN);
                for(double q=q_start;q<q_end+1e-10;q+=q_add){
                    int query_rank = (int) (q * queryN);

                    double merge_sst_v = longToResult(merge_sst_worker.findMinValueWithRank(query_rank));
                    int merge_sst_delta_rank = getDeltaRank(query_a, queryN, merge_sst_v,query_rank);
                    double merge_sst_relative_err = 1.0 * merge_sst_delta_rank / (queryN);
                    err_merge_sst += Math.abs(merge_sst_relative_err) *ratioPerQuery;


                    double full_v = longToResult(full_worker.findMinValueWithRank(query_rank));
                    int full_delta_rank = getDeltaRank(query_a, queryN, full_v, query_rank);
                    double full_relative_err = 1.0 * full_delta_rank / (queryN);
                    err_full += Math.abs(full_relative_err) / (q_count * TEST_CASE);
//                System.out.println("?\t\tfull:"+full_v+" delta:"+full_delta_rank+"\t\tmerge:"+merge_v+" delta:"+merge_delta_rank);
//                    System.out.println("\t\t\trk:"+query_rank+"\t\t"+Math.abs(merge_sst_relative_err));
                }
            }
        }
//        err_merge/=(queryPageNum*pageN);err_full/=(queryPageNum*pageN);
        System.out.println("\t\t" + queryN + "\t" + err_merge_sst + "\t"+err_full);
        //        System.out.println("\t\t\tmerge-point"+"\t"+queryN*(err_mergeBuf-err_full)+"\t"+queryN*(err_merge-err_full));
        err_result.set(RESULT_LINE, err_result.get(RESULT_LINE).concat("\t\t\t\t\t\t" + err_merge_sst + "\t"+err_full));
    }

    public void show_time_result() {
        System.out.println(time_result);
    }

    public static void setTestCase(int tc) {
        TEST_CASE = tc;
    }


    public static void main(String[] args) throws IOException {
        MainForKLLSST main;

        System.out.println("sketch test: byPage & byChunkDivide" + "\n");
        for (int startType=1,endType=2,dataType = startType; dataType <= endType; dataType++){ // CHECK IT
            main = new MainForKLLSST();
            main.prepareA(dataType);
            main.err_result.add("\t");
            for (int queryN : new int[]{40000000/*,chunkN*/})
                for (int query_mem : new int[]{/*1024*16,1024*32,1024*64,*/1024 * 256/*,1024*256*//*,1024*512,1024*1024*/})
                    for(int page_seri : new int[]{4096})
                        for(int SketchSizeRatio:new int[]{1,2,4,8,16})
                        main.testError(queryN, query_mem, page_seri,SketchSizeRatio);
            System.out.println("byPage & byChunkDivide\nTEST_CASE=" + TEST_CASE);
            System.out.println("\nError rate:");
            for (String s : main.err_result)
                System.out.println(s);
        }
    }
}

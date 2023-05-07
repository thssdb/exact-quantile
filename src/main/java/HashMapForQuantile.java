/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * long-long hashmap with <= 65536 elements.
 * rebuild (right shift keys and merge) to keep <=65536.
 */
public class HashMapForQuantile {
  private final int maxSize = (1 << 16) + 1, maxListSize = 1 << 16;
  private final int bucketBits = 16;
  int bitsOfValue, remainingBits, minBits;
  boolean isBucket;long[] bucket;


  int totSize;
  int[] root;
  int[] rc;
  long[] value, count;
  long randSeed;
  private long deltaForUnsignedCompare;
  IntArrayList index;
  long DEBUG;


  public HashMapForQuantile(int bits, int minBits) { // simply 16-bit bucket when remainingBits<=minBits.
    this.minBits = minBits;
    bitsOfValue = remainingBits = bits;
    isBucket = false;

    totSize = 0;
    root = new int[maxSize];
//    rootSize = new long[maxSize];
//    lc = new int[maxSize + 1];  // node id is 1...65536+1
    rc = new int[maxSize + 1];
//    fa = new int[maxSize + 1];
//    rnd = new int[maxSize + 1];
    value = new long[maxSize + 1];
    count = new long[maxSize + 1];
    randSeed = 233L;//System.currentTimeMillis();
    if (bits == 64) deltaForUnsignedCompare = 1L << 63; // unsigned long
    else deltaForUnsignedCompare = 0;
  }

  private long hash(long key) { // Thomas Wang's hash (64bit)
    key = (~key) + (key << 21);
    key = key ^ (key >>> 24);
    key = (key + (key << 3)) + (key << 8);
    key = key ^ (key >>> 14);
    key = (key + (key << 2)) + (key << 4);
    key = key ^ (key >>> 28);
    key = key + (key << 31);
    return key;
  }

  private int nextRandInt() { // xorShift
    randSeed ^= (randSeed << 21);
    randSeed ^= (randSeed >>> 35);
    randSeed ^= (randSeed << 4);
    return (int) randSeed;
  }


  private void turnToBucket() {
//    System.out.println("[turnToBucket]+remaining:"+remainingBits+"  tot:"+totSize);
    isBucket = true;
    if(bucket == null)
      bucket = new long[1 << bucketBits];
    else Arrays.fill(bucket,0);
    for(int i=1;i<=totSize;i++)
      bucket[(int)(value[i]>>> (remainingBits - bucketBits))] += count[i];
    remainingBits = bucketBits;
  }

  /*private void rebuild() { // called when total size == maxSize
    if(index==null){
      index = new IntArrayList(maxSize);
//      index.size(maxSize+1);
      for (int i = 0; i < maxSize; i++) index.add(i+1);
    }

    index.sort((x, y) -> Long.compare(value[x]^deltaForUnsignedCompare, value[y]^deltaForUnsignedCompare));
    long xor;
    int SHR = remainingBits - 16, tmpLen;
    for (int i = 1; i < maxSize; i++) {
      xor = value[index.getInt(i)] ^ value[index.getInt(i-1)];
      if ((xor >>> SHR) > 0) continue;
      tmpLen = 0;
      for (int j = 32; j > 0; j >>= 1)
        if ((xor >>> j) > 0) {
          tmpLen += j;
          xor >>>= j;
        }
      tmpLen += xor;
      SHR = Math.min(SHR, tmpLen);
    }
//    System.out.println("[DEBUG][rebuild] findSHR  :" + SHR);


    if(remainingBits - SHR <=minBits) {
      turnToBucket();
      return;
    }else remainingBits-=SHR;

    Arrays.fill(root, 0);
    Arrays.fill(rc, 0);
    totSize = 0;
    deltaForUnsignedCompare = 0;
//    index.sort((x, y) -> Long.compare(count[y],count[x]));
//    if(lc==null)lc=new int[maxSize+1];
//    for (int i = 1; i <= maxSize; i++)lc[i]=rc[i]=i; // re-use array rc
//    for (int i = 1; i <= maxSize; i++) {
//      int x = index.getInt(i-1), y = rc[x];
//      if (y!=i) {
//        long tmp = value[i];
//        value[i] = value[y];
//        value[y] = tmp;
//        tmp = count[i];
//        count[i] = count[y];
//        count[y] = tmp;
//
//        lc[y] = lc[i];
//        lc[i] = x;
//        rc[lc[y]]=y;
//        rc[x]=i;
//      }
//    }
//    Arrays.fill(lc, 0);
//    System.out.println();
//    for(int i=1;i<=4;i++)System.out.print("  ("+value[i]+","+count[i]+") ");
//    for (int i = 2; i <= maxSize; i++)if(count[i]>count[i-1])DEBUG++;
//    System.out.println("!!! " +DEBUG);
    for (int i = 1; i <= maxSize; i++)
      insertHashMap(value[i] >>> SHR, count[i]);
  }*/

  private void rebuild() { // called when total size == maxSize
//    System.out.println("[rebuild]+remaining:"+remainingBits+"  tot:"+totSize);
    if (remainingBits - 1 <= minBits) {
      turnToBucket();
      return;
    } else remainingBits -= 1;

    Arrays.fill(root, 0);
    Arrays.fill(rc, 0);
    totSize = 0;
    deltaForUnsignedCompare = 0;
    for (int i = 1; i <= maxSize; i++)
      insertHashMap(value[i] >>> 1, count[i]);
  }

  private void insertHashMap(final long num, final long freq) { // as list
    int index = (int) (hash(num) & 65535);
    if (root[index] == 0) {
      root[index] = ++totSize;
      value[totSize] = num;
      count[totSize] = freq;
    } else /*if (rootSize[index] <= 65536) */{ // list when size is small
      int x=root[index], last=0;
      while (value[x] != num && rc[x] != 0) {/*last=x;*/x = rc[x];/*DEBUG++;*/}
      if (value[x] == num) {
        count[x] += freq;
      } else {
        rc[x] = ++totSize;
        value[totSize] = num;
        count[totSize] = freq;
      }/*
      if(last!=0 && count[x]>count[last]){
        long tmp = value[x];value[x]=value[last];value[last]=tmp;
        tmp = count[x];count[x]=count[last];count[last]=tmp;
      }*/
    }
    if (totSize == maxSize)
      rebuild();
  }


  public void insert(long num, long freq) {
    num >>>= bitsOfValue - remainingBits;
    if(isBucket){
      bucket[(int)num]+=freq;
    }else insertHashMap(num, freq);
  }

  public int getRemainingBits() {
    return remainingBits;
  }


  public List<Long> findResultIndex(long K1, long K2) {
    List<Long> result = new ArrayList<>(8);
    long sum = 0;

    if(isBucket){
      for(int i=0;i<(1<<bucketBits);i++){
        sum += bucket[i];
        if (sum >= K1 && result.size() == 0) {
          result.add((long)i);
          result.add(sum - bucket[i]);
        }
        if (sum >= K2 && result.size() == 2) {
          result.add((long)i);
          result.add(sum - bucket[i]);
          break;
        }
      }
    }else {
      if(index==null)index = new IntArrayList(totSize);
      index.size(totSize);
      for(int i=0;i<totSize;i++)index.set(i,i+1);
      index.sort((x, y) -> Long.compare(value[x]^deltaForUnsignedCompare, value[y]^deltaForUnsignedCompare));
      int x;
      for (int i = 0; i < totSize; i++) {
        x = index.getInt(i);
//      System.out.println(count[x] + "  " + value[x]);
        sum += count[x];
        if (sum >= K1 && result.size() == 0) {
          result.add(value[x]);
          result.add(sum - count[x]);
        }
        if (sum >= K2 && result.size() == 2) {
          result.add(value[x]);
          result.add(sum - count[x]);
          break;
        }
      }
    }
    return result;
  }

}

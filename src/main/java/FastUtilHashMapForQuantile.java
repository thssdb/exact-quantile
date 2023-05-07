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
import it.unimi.dsi.fastutil.longs.Long2LongMap;
import it.unimi.dsi.fastutil.longs.Long2LongOpenHashMap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * hppc long-long hashmap with <= 65536 elements.
 * rebuild (right shift keys and merge) to keep <=65536.
 */
public class FastUtilHashMapForQuantile {
  private final int maxSize = (1 << 16)+1;
  private final int bucketBits = 16;
  int bitsOfValue, remainingBits, minBits;
  boolean isBucket;long[] bucket;

  Long2LongOpenHashMap hashMap;
  private long deltaForUnsignedCompare;
  IntArrayList index;
  long DEBUG;


  public FastUtilHashMapForQuantile(int bits, int minBits) { // simply 16-bit bucket when remainingBits<=minBits.
    this.minBits = minBits;
    bitsOfValue = remainingBits = bits;
    isBucket = false;
    if (bits == 64) deltaForUnsignedCompare = 1L << 63; // unsigned long
    else deltaForUnsignedCompare = 0;
    index = new IntArrayList(maxSize);
    hashMap = new Long2LongOpenHashMap(maxSize);
  }

  private void turnToBucket() {
//    System.out.println("[turnToBucket]+remaining:"+remainingBits+"  tot:"+totSize);
    isBucket = true;
    if(bucket == null)
      bucket = new long[1 << bucketBits];
    else Arrays.fill(bucket,0);
//    for (LongLongCursor c = hashMap.cursor(); c.moveNext();) {
//      bucket[(int)(c.key()>>> (remainingBits - bucketBits))] +=c.value();
//    }
    hashMap.forEach((k,v)->bucket[(int)(k>>> (remainingBits - bucketBits))] += v);
    remainingBits = bucketBits;
  }

  private void rebuild() { // called when total size == maxSize
//    System.out.println("[rebuild]+remaining:"+remainingBits+"  tot:"+totSize);
    int SHR = 1;
    if (remainingBits - SHR <= minBits) {
      turnToBucket();
      return;
    }
    deltaForUnsignedCompare = 0;
    Long2LongOpenHashMap newMap = new Long2LongOpenHashMap(maxSize);
    hashMap.forEach((k,v)->newMap.addTo(k>>>1,v));
    while (newMap.size() == maxSize) {
      SHR++;
      if (remainingBits - SHR <= minBits) {
        turnToBucket();
        return;
      }
      newMap.clear();
      final int shr = SHR;
      hashMap.forEach((k,v)->newMap.addTo(k>>>shr,v));
    }
    remainingBits -= SHR;
    hashMap = newMap;
  }

  public void insert(long num, long freq) {
    num >>>= bitsOfValue - remainingBits;
    if(isBucket){
      bucket[(int)num]+=freq;
    }else{
      hashMap.addTo(num,freq);
      if (hashMap.size() == maxSize)
        rebuild();
    }
//    for (LongLongCursor c = hashMap.cursor(); c.moveNext();)
//      System.out.print("("+c.key()+","+c.value()+")");
//    System.out.println();
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
      long[] value, count;
      value = new long[maxSize];
      count = new long[maxSize];
      int totSize = hashMap.size();
      index.size(totSize);
      for(int i=0;i<totSize;i++)index.set(i,i);
      int tmp = 0;
      for(Long2LongMap.Entry entry :hashMap.long2LongEntrySet()) {
        value[tmp]=entry.getLongKey();
        count[tmp]=entry.getLongValue();
        tmp++;
      }
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

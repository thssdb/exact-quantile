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

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
/**
 * This is a bitwise trie (01 trie) with a maximum size.
 */
public class FixedBitwiseTrie {
  private static final int maxLeafSize=1<<16;
  private int trieDepth, maxDepth;
//  int maxSize, cntSize;
  List<List<Node>> nodeList;
  long[] count,value;
  // store all the data in hashmap when
  TreeMap<Long, Long>treeMap;

  private static class Node{
    int l,r;
    public Node(){
      l=r=-1;
    }
  }
  private int getNextIndex(Node node, int depth, int bit){
    List<Node> nextList = nodeList.get(depth+1);
    if(bit==0){
      if(node.l==-1){
        node.l = nextList.size();
        nextList.add(new Node());
      }
      return node.l;
    }else{
      if(node.r==-1){
        node.r = nextList.size();
        nextList.add(new Node());
      }
      return node.r;
    }
  }

  // maxDepth = 32 or 64
  public FixedBitwiseTrie(int maxDepth){
    trieDepth = this.maxDepth = maxDepth;
    treeMap = new TreeMap<>();
  }
  private void buildTrie(){
    count = new long[maxLeafSize];
    value = new long[maxLeafSize];
    nodeList = new ArrayList<>(trieDepth+1);
    for(int i=0;i<=trieDepth;i++)
      nodeList.add(new ArrayList<>(/*Math.min(1<<i, maxLeafSize)*/));
    nodeList.get(0).add(new Node()); // root
//    System.out.println("   mid mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
    for(Map.Entry<Long, Long> entry:treeMap.entrySet()){
      insertToTrie(entry.getKey(),entry.getValue());
    }
//    System.out.println("   mid mem:"+(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024/1024);
  }
  private void mergeLeaf(){
    nodeList.remove(trieDepth);
    List<Node> faList = nodeList.get(trieDepth-1);
    Node faNode;
    for(int i=0;i<faList.size();i++) {
      faNode = faList.get(i);
      if(faNode.l!=-1&&faNode.r!=-1) {
        count[i] = count[faNode.l] + count[faNode.r];
        value[i] = value[faNode.l]>>>1;
      }else{
        int nextIndex = faNode.l+faNode.r+1;
        count[i]=count[nextIndex];
        value[i]=value[nextIndex]>>>1;
      }
    }
    for(int i=faList.size();i<maxLeafSize;i++)
      count[i]=value[i]=0;
  }
  private void insertToTrie(long index, long times){
    int i, bit, nextIndex=-1;
    Node cntNode = nodeList.get(0).get(0);
    for(i=1;i<=trieDepth;i++){
      bit = (int)((index>>>(trieDepth-i))&1);
//      cntIndex = getNextIndex(nodeList.get(i-1).get(cntIndex), i-1, bit);
//      int nextIndex = bit==0?nodeList.get(i-1).get(cntIndex).l:nodeList.get(i-1).get(cntIndex).r;
      nextIndex = bit==0?cntNode.l:cntNode.r;
      if(nextIndex == -1)break;
      else cntNode=nodeList.get(i).get(nextIndex);
    }
    if(i>trieDepth){
      // exist in trie
      count[nextIndex]+=times;
      return;
    }
    while(nodeList.get(trieDepth).size()==maxLeafSize) {
      mergeLeaf();
      trieDepth-=1;
      index>>>=1;
    }
    for(;i<=trieDepth;i++){
      bit = (int)((index>>>(trieDepth-i))&1);
      if(bit==0)
        cntNode.l=nodeList.get(i).size();
      else cntNode.r=nodeList.get(i).size();
      cntNode = new Node();
      nodeList.get(i).add(cntNode);
    }
    nextIndex = nodeList.get(trieDepth).size()-1;
    count[nextIndex] = times;
    value[nextIndex] = index;
  }

  public void insert(long maskedLongBits, long times){
    long index = maskedLongBits>>>(maxDepth-trieDepth);
    if(treeMap!=null){
      treeMap.merge(index, times, Long::sum);
      if(treeMap.size()>=maxLeafSize){
        buildTrie();
        treeMap = null;
      }
    }else insertToTrie(index, times);
  }
  public List<Long> findResultIndex(long K1,long K2){
    List<Long> result = new ArrayList<>(8);
    long sum = 0;
    if(treeMap!=null){
      for(Map.Entry<Long, Long> entry:treeMap.entrySet()){
        sum+=entry.getValue();
        if(sum>=K1&&result.size()==0){
          result.add(entry.getKey());
          result.add(sum-entry.getValue());
        }
        if(sum>=K2&&result.size()==2){
          result.add(entry.getKey());
          result.add(sum-entry.getValue());
          break;
        }
      }
    }else {
      for (int i = 0; i < nodeList.get(trieDepth).size(); i++) {
//        System.out.println(count[i]+"  "+value[i]);
        sum += count[i];
        if(sum>=K1&&result.size()==0){
          result.add(value[i]);
          result.add(sum-count[i]);
        }
        if(sum>=K2&&result.size()==2){
          result.add(value[i]);
          result.add(sum-count[i]);
          break;
        }
      }
    }
    return result;
  }
  public int getTrieDepth(){return trieDepth;}
}

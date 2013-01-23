/* 
 * File:   outputsorter.h
 * Author: enckevort
 *
 * Created on 30 mei 2011, 13:24
 */

#ifndef OUTPUTSORTER_H
#define	OUTPUTSORTER_H

#include <vector>

class OutputSorter
{
public:
  OutputSorter (const unsigned &NumBoxes, const std::string &CurrentChr,
                std::ofstream &InvOutf);
  void SortAndOutputInversions (std::vector<SPLIT_READ> &Reads, std::vector<unsigned> Inv[]);
  void SortAndOutputNonTemplateInversions (std::vector<SPLIT_READ> &Reads, std::vector<unsigned> Inv[]);
  virtual ~OutputSorter ();
private:
  OutputSorter (const OutputSorter&);
  int ReportIndelEvents (std::vector<Indel4output> &IndelEvents,
                         std::vector<SPLIT_READ> &GoodIndels);
  
  int DoSortAndOutputInversions (std::vector<SPLIT_READ> &Reads,
                                 std::vector<unsigned> Inv[],
                                 bool isNonTemplateInversion);
  unsigned NumBoxes;
  std::string* CurrentChr;
  std::ofstream* InvOutf;
};

#endif	/* OUTPUTSORTER_H */


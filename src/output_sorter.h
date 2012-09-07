/* 
 * File:   outputsorter.h
 * Author: enckevort
 *
 * Created on 30 mei 2011, 13:24
 */

#ifndef OUTPUTSORTER_H
#define	OUTPUTSORTER_H

#include <vector>
#include "control_state.h"

class OutputSorter
{
public:
  OutputSorter (const unsigned &NumBoxes, const std::string &CurrentChr,
                std::ofstream &InvOutf);
  void SortAndOutputInversions (ControlState& currentState, std::vector<SPLIT_READ> &Reads, std::vector<unsigned> Inv[]);
  void SortAndOutputNonTemplateInversions (ControlState& currentState, std::vector<SPLIT_READ> &Reads, std::vector<unsigned> Inv[]);
  virtual ~OutputSorter ();
private:
  OutputSorter (const OutputSorter&);
  int ReportIndelEvents (ControlState& currentState,
                         std::vector<Indel4output> &IndelEvents,
                         std::vector<SPLIT_READ> &GoodIndels);
  
  int DoSortAndOutputInversions (ControlState& currentState,
                                 std::vector<SPLIT_READ> &Reads,
                                 std::vector<unsigned> Inv[],
                                 bool isNonTemplateInversion);
  unsigned NumBoxes;
  std::string* CurrentChr;
  std::ofstream* InvOutf;
};

#endif	/* OUTPUTSORTER_H */


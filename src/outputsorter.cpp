/* 
 * File:   outputsorter.cpp
 * Author: enckevort
 * 
 * Created on 30 mei 2011, 13:24
 */
#include <iostream>
#include <fstream>

#include "pindel.h"
#include "logdef.h"
#include "outputsorter.h"

OutputSorter::OutputSorter (const unsigned &NumBoxes_in, const std::string & CurrentChr_in,
							std::ofstream & InvOutf_in)
{
  NumBoxes = NumBoxes_in;
  CurrentChr = new std::string (CurrentChr_in);
  InvOutf = &InvOutf_in;
}

OutputSorter::~OutputSorter () { }

void
OutputSorter::SortAndOutputInversions (std::vector < SPLIT_READ > &Reads, std::vector<unsigned> Inv[])
{
  LOG_INFO (std::cout << "Sorting and outputing Inversions ..." << std::endl);
  DoSortAndOutputInversions(Reads, Inv, false);
  LOG_INFO (std::cout << "Inversions (INV): " << NumberOfInvInstances << std::endl << std::endl);
}

void
OutputSorter::SortAndOutputNonTemplateInversions (std::vector<SPLIT_READ>& Reads, std::vector<unsigned> Inv[])
{
  LOG_INFO (std::cout << "Sorting and outputing Inversions with non-template sequence ..." <<
			std::endl);
  int Count_INV_NT_output = DoSortAndOutputInversions(Reads, Inv, true);

  LOG_INFO (std::cout << "Inversions with non-template sequence (INV_NT): " <<
			Count_INV_NT_output << std::endl << std::endl);
}

int
OutputSorter::DoSortAndOutputInversions (std::vector<SPLIT_READ> &Reads,
										 std::vector<unsigned> Inv[],
										 bool isNonTemplateInversion)
{
  unsigned int InversionsNum;
  short CompareResult;
  unsigned Temp4Exchange;
  unsigned int GoodNum;
  std::vector<SPLIT_READ> GoodIndels;
  std::vector<Indel4output> IndelEvents;
  int ReportedEventCount = 0;

  for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
	{
	  if (Inv[Box_index].size () >= NumRead2ReportCutOff)
		{
		  InversionsNum = Inv[Box_index].size ();
		  LOG_DEBUG (std::cout << Box_index << "\t" << Inv[Box_index].size () << std::endl);
		  for (unsigned int First = 0; First < InversionsNum - 1; First++)
			{
			  for (unsigned int Second = First + 1; Second < InversionsNum;
				  Second++)
				{
				  LOG_DEBUG (std::cout << InputIndels[First].BPLeft << "\t" << InputIndels[First].BPRight << "\t"
							 InputIndels[Second].BPLeft << "\t" << InputIndels[Second].BPRight << std::endl);
				  {
					if (Reads[Inv[Box_index][First]].ReadLength ==
						Reads[Inv[Box_index][Second]].ReadLength)
					  {
						if (Reads[Inv[Box_index][First]].LeftMostPos ==
							Reads[Inv[Box_index][Second]].LeftMostPos)
						  Reads[Inv[Box_index][Second]].Unique = false;
					  }
					if (Reads[Inv[Box_index][First]].BPLeft <
						Reads[Inv[Box_index][Second]].BPLeft)
					  continue;
					else if (Reads[Inv[Box_index][First]].BPLeft >
							 Reads[Inv[Box_index][Second]].BPLeft)
					  {
						CompareResult = 1;
					  }
					else if (Reads[Inv[Box_index][First]].BPLeft ==
							 Reads[Inv[Box_index][Second]].BPLeft)
					  {
						if (Reads[Inv[Box_index][First]].BPRight <
							Reads[Inv[Box_index][Second]].BPRight)
						  continue;
						else if (Reads[Inv[Box_index][First]].BPRight >
								 Reads[Inv[Box_index][Second]].BPRight)
						  {
							CompareResult = 1;
						  }
						else if (isNonTemplateInversion)
						  {
							if (Reads[Inv[Box_index][First]].NT_size <
								Reads[Inv[Box_index][Second]].NT_size)
							  continue;
							else if (Reads[Inv[Box_index][First]].NT_size >
									 Reads[Inv[Box_index][Second]].NT_size)
							  CompareResult = 1;
						  }
					  }
					if (CompareResult == 1)
					  {
						Temp4Exchange = Inv[Box_index][First];
						Inv[Box_index][First] = Inv[Box_index][Second];
						Inv[Box_index][Second] = Temp4Exchange;
					  }
				  }
				}
			}
		  GoodIndels.clear ();
		  IndelEvents.clear ();

		  for (unsigned int First = 0; First < InversionsNum; First++)
			{
			  GoodIndels.push_back (Reads[Inv[Box_index][First]]);
			}

		  GoodNum = GoodIndels.size ();
		  LOG_DEBUG (std::cout << "GoodNum " << Box_index << " " << GoodNum << std::endl);
		  if (GoodNum == 0)
			continue;
		  LOG_DEBUG (std::cout << GoodNum << std::endl);
		  Indel4output OneIndelEvent;
		  OneIndelEvent.Start = 0;
		  OneIndelEvent.End = 0;
		  OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
		  OneIndelEvent.BPLeft = GoodIndels[0].BPLeft;
		  OneIndelEvent.BPRight = GoodIndels[0].BPRight;
		  OneIndelEvent.WhetherReport = true;
		  for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++)
			{
			  LOG_DEBUG (std::cout << GoodIndex <<
					  "\t" << GoodIndels[GoodIndex].BPLeft <<
					  "\t" << GoodIndels[GoodIndex].BPRight <<
					  "\t" << OneIndelEvent.BPLeft <<
					  "\t" << OneIndelEvent.BPRight << std::endl);
			  if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
				  && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
				OneIndelEvent.End = GoodIndex;
			  else
				{
				  OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
				  OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
				  OneIndelEvent.Support =
					  OneIndelEvent.End - OneIndelEvent.Start + 1;
				  IndelEvents.push_back (OneIndelEvent);
				  OneIndelEvent.Start = GoodIndex;
				  OneIndelEvent.End = GoodIndex;
				  OneIndelEvent.BPLeft = GoodIndels[GoodIndex].BPLeft;
				  OneIndelEvent.BPRight = GoodIndels[GoodIndex].BPRight;
				}
			}

		  OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
		  OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
		  OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
		  LOG_DEBUG (std::cout << OneIndelEvent.Support
				  << "\t" << OneIndelEvent.Start
				  << "\t" << OneIndelEvent.End << std::endl);
		  IndelEvents.push_back (OneIndelEvent);
		  LOG_DEBUG (std::cout << "IndelEvent: " << IndelEvents.size () << std::endl);

		  if (IndelEvents.size ())
			{
			  ReportedEventCount += ReportIndelEvents (IndelEvents, GoodIndels);
			}
		}
	}
  return ReportedEventCount;
}

int
OutputSorter::ReportIndelEvents (std::vector<Indel4output> &IndelEvents,
								 std::vector<SPLIT_READ> &GoodIndels)
{
  int ReportedEventCount = 0;
  for (unsigned EventIndex = 0; EventIndex < IndelEvents.size ();
	  EventIndex++)
	{
	  LOG_DEBUG (std::cout << IndelEvents[EventIndex].Start <<
			  "\t" << IndelEvents[EventIndex].End <<
			  "\t" << IndelEvents[EventIndex].Support << std::endl);
	  unsigned int RealStart = IndelEvents[EventIndex].RealStart;
	  unsigned int RealEnd = IndelEvents[EventIndex].RealEnd;
	  if (IndelEvents[EventIndex].Support < NumRead2ReportCutOff)
		continue;
	  if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize <
		  BalanceCutoff)
		{
		  OutputInversions (GoodIndels, *CurrentChr,
						  IndelEvents[EventIndex].Start,
						  IndelEvents[EventIndex].End,
						  RealStart, RealEnd, *InvOutf);
		  NumberOfInvInstances++;
		  ReportedEventCount++;
		}
	  else
		if (ReportEvent
			(GoodIndels, IndelEvents[EventIndex].Start,
			IndelEvents[EventIndex].End))
		{
		  OutputInversions (GoodIndels, *CurrentChr,
						  IndelEvents[EventIndex].Start,
						  IndelEvents[EventIndex].End,
						  RealStart, RealEnd, *InvOutf);
		  NumberOfInvInstances++;
		  ReportedEventCount++;
		}
	}
  return ReportedEventCount;
}

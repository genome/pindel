#include <iostream>
#include <fstream>
#include <cmath>
// Pindel
#include "pindel.h"
#include "reader.h"
#include "searcher.h"
#include "reporter.h"

using namespace std;

void OutputTDs(const vector <SPLIT_READ> & TDs,
					const string & TheInput, 
					const unsigned int & C_S,
					const unsigned int & C_E,
					const unsigned int & RealStart,
					const unsigned int & RealEnd,
					ofstream & TDOutf) {
   //short ReadLength = Deletions[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;
	
   SupportPerSample NumSupportPerTag[VectorTag.size()];
	
   for (short i = 0; i < VectorTag.size(); i++) {
		NumSupportPerTag[i].NumPlus = 0;
		NumSupportPerTag[i].NumMinus = 0;
		NumSupportPerTag[i].NumUPlus = 0;
		NumSupportPerTag[i].NumUMinus = 0;
	}
	
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	      if (TDs[i].Tag == VectorTag[j]) {
				if (TDs[i].MatchedD == Plus) {
					NumSupportPerTag[j].NumPlus++;
					if (TDs[i].Unique)
					   NumSupportPerTag[j].NumUPlus++;
				}
				else {
					NumSupportPerTag[j].NumMinus++;
					if (TDs[i].Unique)
					   NumSupportPerTag[j].NumUMinus++;
				}
			}
      }
      if (TDs[i].MatchedD == Plus) {
         LeftScore += TDs[i].score;
         LeftS++;
			if (TDs[i].Unique)
			   LeftUNum++;
      }
      else {
         RightScore += TDs[i].score;
         RightS++;
			if (TDs[i].Unique)
			   RightUNum++;
      }
   }
	
	short NumberSupSamples = 0; 
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
   for (short i = 0; i < VectorTag.size(); ++i) {
		if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) ++NumberSupSamples;
		if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) ++NumU_SupSamples;
		Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
	}
	//  << "\t" << Num_U_Reads
	//  << "\t" << NumU_SupSamples
	
   unsigned int EasyScore = LeftS * RightS;
   double PreciseScore = (LeftScore + RightScore) * (-1);
   short GapSize =0;
   if (TDs[C_S].IndelSize < 14) GapSize = TDs[C_S].IndelSize;
   else GapSize = 13 + (int)log10(TDs[C_S].IndelSize - 10);
	CurrentChrMask[TDs[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[TDs[C_S].BPRight + SpacerBeforeAfter] = 'B';
	
   TDOutf << "####################################################################################################" << endl;
   TDOutf << NumberOfTDInstances << "\tTD " << TDs[C_S].IndelSize// << " bases " 
	<< "\tNT " << TDs[C_S].NT_size << " \"" << TDs[C_S].NT_str << "\""
	<< "\tChrID " << TDs[C_S].FragName 
	<< "\tBP " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2
	<< "\tBP_range " << TDs[C_S].BPLeft << "\t" << TDs[C_S].BPRight + 2
	<< "\tSupports " << NumberOfReads  << "\t" << Num_U_Reads
	<< "\t+ " << LeftS - 1  << "\t" << LeftUNum
	<< "\t- " << RightS - 1  << "\t" << RightUNum
	<< "\tS1 " << EasyScore ; //EWL070111 << "\tS2 " << PreciseScore; 
	
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += TDs[i].MS;
   TDOutf << "\tSUM_MS " << SUM_MS;
	
   TDOutf << "\t" << VectorTag.size() << "\tNumSupSamples " << NumberSupSamples << "\t" << NumU_SupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		TDOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
		<< " " << NumSupportPerTag[i].NumUPlus
		<< " " << NumSupportPerTag[i].NumMinus
		<< " " << NumSupportPerTag[i].NumUMinus;
   TDOutf << endl;
	
   //DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength                 
   TDOutf << TheInput.substr(TDs[C_S].BPRight + SpacerBeforeAfter - ReportLength + 1, ReportLength);// << endl;// 
	if (TDs[C_S].NT_size) {
		for (short i = 0; i < TDs[C_S].NT_size; i++) TDOutf << " ";
	}
	//ReportLength 
   TDOutf << Cap2Low(TheInput.substr(TDs[C_S].BPLeft + SpacerBeforeAfter, ReportLength)) << endl;
	
   short SpaceBeforeReadSeq;    
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) { 
      SpaceBeforeReadSeq = ReportLength - TDs[GoodIndex].BP - 1;
      
      for (int i = 0; i < SpaceBeforeReadSeq; i++) TDOutf << " ";
      short SpaceBeforeD = ReportLength + ReportLength - SpacerBeforeAfter - TDs[GoodIndex].ReadLength;
      if (TDs[GoodIndex].MatchedD == Minus) {
         TDOutf << TDs[GoodIndex].UnmatchedSeq << endl;
         //for (int i = 0; i < GapSize; i++) TDOutf << " ";    
         //TDOutf << TDs[GoodIndex].UnmatchedSeq.substr(TDs[GoodIndex].BP + 1, TDs[GoodIndex].ReadLength - TDs[GoodIndex].BP);// << endl;
      }
      else {
         TDOutf << ReverseComplement(TDs[GoodIndex].UnmatchedSeq) << endl; 
         //for (int i = 0; i < GapSize; i++) TDOutf << " ";  
         //TDOutf << ReverseComplement(TDs[GoodIndex].UnmatchedSeq).substr(TDs[GoodIndex].BP + 1, TDs[GoodIndex].ReadLength - TDs[GoodIndex].BP);// << endl;
      } 
      //for (int i = 0; i < SpaceBeforeD; i++) TDOutf << " ";
      TDOutf << "\t" << TDs[GoodIndex].MatchedD << "\t" 
		<< TDs[GoodIndex].MatchedRelPos 
		<< "\t" << TDs[GoodIndex].MS
		<< "\t" << TDs[GoodIndex].Tag 
		<< "\t" <<  TDs[GoodIndex].Name << endl;
		//<< "\t" << TDs[C_S].BPLeft
		//<< "\t" << TDs[C_S].BPRight << endl; 
   }     
}

void OutputDeletions(const vector <SPLIT_READ> & Deletions,
                     const string & TheInput, 
                     const unsigned int & C_S,
                     const unsigned int & C_E,
							const unsigned int & RealStart,
							const unsigned int & RealEnd,
                     ofstream & DeletionOutf) {
   //short ReadLength = Deletions[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;
	//cout << "d_1" << endl;
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;
	//cout << "d_2" << endl;
   SupportPerSample NumSupportPerTag[VectorTag.size()];
	
   for (short i = 0; i < VectorTag.size(); i++) {
		NumSupportPerTag[i].NumPlus = 0;
		NumSupportPerTag[i].NumMinus = 0;
		NumSupportPerTag[i].NumUPlus = 0;
		NumSupportPerTag[i].NumUMinus = 0;
	}
	//cout << "d_3" << endl;
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	      if (Deletions[i].Tag == VectorTag[j]) {
				if (Deletions[i].MatchedD == Plus) {
					NumSupportPerTag[j].NumPlus++;
					if (Deletions[i].Unique)
					   NumSupportPerTag[j].NumUPlus++;
				}
				else {
					NumSupportPerTag[j].NumMinus++;
					if (Deletions[i].Unique)
					   NumSupportPerTag[j].NumUMinus++;
				}
			}
      }
      if (Deletions[i].MatchedD == Plus) {
         LeftScore += Deletions[i].score;
         LeftS++;
			if (Deletions[i].Unique)
			   LeftUNum++;
      }
      else {
         RightScore += Deletions[i].score;
         RightS++;
			if (Deletions[i].Unique)
			   RightUNum++;
      }
   }
	//cout << "d_4" << endl;
	short NumberSupSamples = 0; 
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
   for (short i = 0; i < VectorTag.size(); ++i) {
		if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) ++NumberSupSamples;
		if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) ++NumU_SupSamples;
		Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
	}
	//  << "\t" << Num_U_Reads
	//cout << "d_5" << endl;
   unsigned int EasyScore = LeftS * RightS;
   double PreciseScore = (LeftScore + RightScore) * (-1);
   short GapSize = 0;
   if (Deletions[C_S].IndelSize < 14) GapSize = Deletions[C_S].IndelSize;
   else GapSize = 13 + (int)log10(Deletions[C_S].IndelSize - 10);
	//cout << "d_5a" << endl;
	CurrentChrMask[Deletions[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	//cout << "d_5b" << endl;
	CurrentChrMask[Deletions[C_S].BPRight + SpacerBeforeAfter] = 'B';
	//cout << "d_5c" << endl;
	CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B';
	//cout << "d_5d" << endl;
	CurrentChrMask[RealEnd + SpacerBeforeAfter] = 'B';
	//cout << "d_5e" << endl;
   DeletionOutf << "####################################################################################################" << endl;
   DeletionOutf << NumberOfDeletionsInstances << "\tD " << Deletions[C_S].IndelSize// << " bases " 
	<< "\tNT " << Deletions[C_S].NT_size << " \"" << Deletions[C_S].NT_str << "\""
	<< "\tChrID " << Deletions[C_S].FragName 
	<< "\tBP " << Deletions[C_S].BPLeft + 1 << "\t" << Deletions[C_S].BPRight + 1 
	<< "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1  
	<< "\tSupports " << NumberOfReads << "\t" << Num_U_Reads
	<< "\t+ " << LeftS - 1 << "\t" << LeftUNum
	<< "\t- " << RightS - 1 << "\t" << RightUNum 
	<< "\tS1 " << EasyScore; //EWL070111  << "\tS2 " << PreciseScore; 
	//cout << "d_6" << endl;
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += Deletions[i].MS;
   DeletionOutf << "\tSUM_MS " << SUM_MS;
	
   DeletionOutf << "\t" << VectorTag.size() << "\tNumSupSamples " << NumberSupSamples << "\t" << NumU_SupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		DeletionOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
		<< " " << NumSupportPerTag[i].NumUPlus
		<< " " << NumSupportPerTag[i].NumMinus
		<< " " << NumSupportPerTag[i].NumUMinus;
   DeletionOutf << endl;
	//cout << "d_7" << endl;
   //DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength                 
   DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength);// << endl;// ReportLength    
   if (Deletions[C_S].IndelSize >= 14) {
      DeletionOutf << Cap2Low(TheInput.substr(Deletions[C_S].Left + Deletions[C_S].BP + 1, 5))
		<< "<" << Deletions[C_S].IndelSize - 10 << ">" 
		<< Cap2Low(TheInput.substr(Deletions[C_S].Right - Deletions[C_S].ReadLength + Deletions[C_S].BP - 3, 5));  
   }   
   else DeletionOutf << Cap2Low(TheInput.substr(Deletions[C_S].Left + Deletions[C_S].BP + 1, GapSize));             
   DeletionOutf << TheInput.substr(Deletions[C_S].Left + Deletions[C_S].BP + 1 + Deletions[C_S].IndelSize, ReportLength - GapSize) << endl;// ReportLength
   short SpaceBeforeReadSeq;    
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) { 
      SpaceBeforeReadSeq = ReportLength - Deletions[GoodIndex].BP - 1;
      
      for (int i = 0; i < SpaceBeforeReadSeq; i++) DeletionOutf << " ";
      short SpaceBeforeD = ReportLength + ReportLength - SpaceBeforeReadSeq - Deletions[GoodIndex].ReadLength;
      if (Deletions[GoodIndex].MatchedD == Minus) {
         DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr(0, Deletions[GoodIndex].BP + 1);// << endl;
         for (int i = 0; i < GapSize; i++) DeletionOutf << " ";    
         DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      }
      else {
         DeletionOutf << ReverseComplement(Deletions[GoodIndex].UnmatchedSeq).substr(0, Deletions[GoodIndex].BP + 1);// << endl; 
         for (int i = 0; i < GapSize; i++) DeletionOutf << " ";  
         DeletionOutf << ReverseComplement(Deletions[GoodIndex].UnmatchedSeq).substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      } 
      for (int i = 0; i < SpaceBeforeD; i++) DeletionOutf << " ";
      DeletionOutf << "\t" << Deletions[GoodIndex].MatchedD << "\t" 
		<< Deletions[GoodIndex].MatchedRelPos 
		<< "\t" << Deletions[GoodIndex].MS
		<< "\t" << Deletions[GoodIndex].Tag 
		<< "\t" <<  Deletions[GoodIndex].Name << endl;
		//<< "\t" << Deletions[C_S].BPLeft
		//<< "\t" << Deletions[C_S].BPRight << endl; 
   }     
}

void OutputInversions(const vector <SPLIT_READ> & Inv, 
							 const string & TheInput, 
							 const unsigned int & C_S, 
							 const unsigned int & C_E,
							 const unsigned int & RealStart,
							 const unsigned int & RealEnd,
							 ofstream & InvOutf) {
	int LeftNT_index = -1;
	int RightNT_index = -1;
	for (unsigned Index = C_S; Index <= C_E; Index++) {
		if (Inv[Index].MatchedD == Plus) {
			LeftNT_index = Index;
			break;
		}
	}
	for (unsigned Index = C_S; Index <= C_E; Index++) {
		if (Inv[Index].MatchedD == Minus) {
			RightNT_index = Index;
			break;
		}
	}
	short LeftNT_size = 0;
	short RightNT_size = 0;
	string LeftNT_str = "";
	string RightNT_str = "";
	if (LeftNT_index != -1) {
		LeftNT_size = Inv[LeftNT_index].NT_size;
		LeftNT_str = Inv[LeftNT_index].NT_str;
	}
	if (RightNT_index != -1) {
		RightNT_size = Inv[RightNT_index].NT_size;
		RightNT_str = Inv[RightNT_index].NT_str;
	}
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;
	
   SupportPerSample NumSupportPerTag[VectorTag.size()];
	
   for (short i = 0; i < VectorTag.size(); i++) {
		NumSupportPerTag[i].NumPlus = 0;
		NumSupportPerTag[i].NumMinus = 0;
		NumSupportPerTag[i].NumUPlus = 0;
		NumSupportPerTag[i].NumUMinus = 0;
	}
	
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	      if (Inv[i].Tag == VectorTag[j]) {
				if (Inv[i].MatchedD == Plus) {
					NumSupportPerTag[j].NumPlus++;
					if (Inv[i].Unique)
					   NumSupportPerTag[j].NumUPlus++;
				}
				else {
					NumSupportPerTag[j].NumMinus++;
					if (Inv[i].Unique)
					   NumSupportPerTag[j].NumUMinus++;
				}
			}
      }
      if (Inv[i].MatchedD == Plus) {
         LeftScore += Inv[i].score;
         LeftS++;
			if (Inv[i].Unique)
			   LeftUNum++;
      }
      else {
         RightScore += Inv[i].score;
         RightS++;
			if (Inv[i].Unique)
			   RightUNum++;
      }
   }
	
	short NumberSupSamples = 0; 
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
   for (short i = 0; i < VectorTag.size(); ++i) {
		if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) ++NumberSupSamples;
		if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) ++NumU_SupSamples;
		Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
	}
	//  << "\t" << Num_U_Reads
	//  << "\t" << NumU_SupSamples
	
   unsigned int EasyScore = LeftS * RightS;
   //double PreciseScore = (LeftScore + RightScore) * (-1);
   short GapSize =0;
   if (Inv[C_S].IndelSize < 14) GapSize = Inv[C_S].IndelSize;
   else GapSize = 13 + (int)log10(Inv[C_S].IndelSize - 10);
	CurrentChrMask[Inv[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[Inv[C_S].BPRight + SpacerBeforeAfter] = 'B';
   InvOutf << "####################################################################################################" << endl;
   InvOutf << NumberOfInvInstances << "\tINV " << Inv[C_S].IndelSize// << " bases " 
	<< "\tNT " <<  LeftNT_size << ":" << RightNT_size << " \"" << LeftNT_str << "\":\"" << RightNT_str << "\""
	<< "\tChrID " << Inv[C_S].FragName 
	<< "\tBP " << Inv[C_S].BPLeft + 1 - 1 << "\t" << Inv[C_S].BPRight + 1 + 1 
	<< "\tBP_range " << Inv[C_S].BPLeft + 1 - 1 << "\t" << Inv[C_S].BPRight + 1 + 1
	<< "\tSupports " << NumberOfReads  << "\t" << Num_U_Reads
	<< "\t+ " << LeftS - 1  << "\t" << LeftUNum
	<< "\t- " << RightS - 1  << "\t" << RightUNum
	<< "\tS1 " << EasyScore; //EWL070111  << "\tS2 " << PreciseScore; 
	
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += Inv[i].MS;
   InvOutf << "\tSUM_MS " << SUM_MS;
	
   InvOutf << "\t" << VectorTag.size() << "\tNumSupSamples " << NumberSupSamples << "\t" << NumU_SupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		InvOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
		<< " " << NumSupportPerTag[i].NumUPlus
		<< " " << NumSupportPerTag[i].NumMinus
		<< " " << NumSupportPerTag[i].NumUMinus;
   InvOutf << endl;
	
	short SpaceBeforeReadSeq;  
   //DeletionOutf << TheInput.substr(Deletions[C_S].Left - ReportLength + Deletions[C_S].BP + 1, ReportLength * 2) << endl;// << endl;// ReportLength                 
   InvOutf << TheInput.substr(Inv[C_S].BPLeft + SpacerBeforeAfter - ReportLength, ReportLength);//;// << endl;// ReportLength  
	//cout << Inv[C_S].NT_size << "\t" << Inv[C_S].NT_2size << endl;
	if (LeftNT_size) {
		for (int i = 0; i < LeftNT_size; i++) {
			InvOutf << " ";
		}
	}
	InvOutf << Cap2Low( ReverseComplement(TheInput.substr(Inv[C_S].BPRight + 1 + SpacerBeforeAfter - ReportLength, ReportLength))) << endl; 
	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) { 
      if (Inv[GoodIndex].MatchedD == Plus && Inv[GoodIndex].MatchedRelPos < Inv[GoodIndex].BPLeft) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << ReverseComplement(Inv[GoodIndex].UnmatchedSeq);
			for (int i = 0; i < Inv[GoodIndex].BP; i++) InvOutf << " ";
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
			<< "\t" << Inv[GoodIndex].MatchedD << "\t" 
			<< Inv[GoodIndex].MatchedRelPos 
			<< "\t" << Inv[GoodIndex].MS
			<< "\t" << Inv[GoodIndex].Tag 
			<< "\t" <<  Inv[GoodIndex].Name << endl;
			//<< "\t" << Deletions[C_S].BPLeft
			//<< "\t" << Deletions[C_S].BPRight << endl; 
      }
      else if (Inv[GoodIndex].MatchedD == Plus && Inv[GoodIndex].MatchedRelPos > Inv[GoodIndex].BPLeft) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << Inv[GoodIndex].UnmatchedSeq;
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
			<< "\t" << Inv[GoodIndex].MatchedD << "\t" 
			<< Inv[GoodIndex].MatchedRelPos 
			<< "\t" << Inv[GoodIndex].MS
			<< "\t" << Inv[GoodIndex].Tag 
			<< "\t" <<  Inv[GoodIndex].Name << endl;
      } 
   }  
	
	InvOutf << "----------------------------------------------------------------------------------------------------" << endl;
	
	
	InvOutf << Cap2Low( ReverseComplement(TheInput.substr(Inv[C_S].BPLeft + SpacerBeforeAfter, ReportLength)));
	if (RightNT_size) {
		for (int i = 0; i < RightNT_size; i++) {
			InvOutf << " ";
		}
	}
	InvOutf << TheInput.substr(Inv[C_S].BPRight + 1 + SpacerBeforeAfter, ReportLength) << endl; 
	for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) { 
      if (Inv[GoodIndex].MatchedD == Minus && Inv[GoodIndex].MatchedRelPos < Inv[GoodIndex].BPRight) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << ReverseComplement(Inv[GoodIndex].UnmatchedSeq);
			for (int i = 0; i < Inv[GoodIndex].BP; i++) InvOutf << " ";
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str 
			<< "\t" << Inv[GoodIndex].MatchedD << "\t" 
			<< Inv[GoodIndex].MatchedRelPos 
			<< "\t" << Inv[GoodIndex].MS
			<< "\t" << Inv[GoodIndex].Tag 
			<< "\t" <<  Inv[GoodIndex].Name << endl;
			//<< "\t" << Deletions[C_S].BPLeft
			//<< "\t" << Deletions[C_S].BPRight << endl; 
      }
      else if (Inv[GoodIndex].MatchedD == Minus && Inv[GoodIndex].MatchedRelPos > Inv[GoodIndex].BPRight) {
			SpaceBeforeReadSeq = ReportLength - Inv[GoodIndex].BP - 1;
			for (int i = 0; i < SpaceBeforeReadSeq; i++) InvOutf << " ";
         InvOutf << Inv[GoodIndex].UnmatchedSeq;
			InvOutf //<< "\t" << Inv[GoodIndex].NT_size << "\t\"" << Inv[GoodIndex].NT_str
			<< "\t" << Inv[GoodIndex].MatchedD << "\t" 
			<< Inv[GoodIndex].MatchedRelPos 
			<< "\t" << Inv[GoodIndex].MS
			<< "\t" << Inv[GoodIndex].Tag 
			<< "\t" <<  Inv[GoodIndex].Name << endl;
      } 
   } 
	
	
}

void OutputSIs(const vector <SPLIT_READ> & SIs,
               const string & TheInput,
               const unsigned int & C_S,
               const unsigned int & C_E,
					const unsigned int & RealStart,
					const unsigned int & RealEnd,
               ofstream & SIsOutf) {
   //short ReadLength = SIs[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;               
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;
	
   SupportPerSample NumSupportPerTag[VectorTag.size()];
	
   for (short i = 0; i < VectorTag.size(); i++) {
		NumSupportPerTag[i].NumPlus = 0;
		NumSupportPerTag[i].NumMinus = 0;
		NumSupportPerTag[i].NumUPlus = 0;
		NumSupportPerTag[i].NumUMinus = 0;
	}
	
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	      if (SIs[i].Tag == VectorTag[j]) {
				if (SIs[i].MatchedD == Plus) {
					NumSupportPerTag[j].NumPlus++;
					if (SIs[i].Unique)
					   NumSupportPerTag[j].NumUPlus++;
				}
				else {
					NumSupportPerTag[j].NumMinus++;
					if (SIs[i].Unique)
					   NumSupportPerTag[j].NumUMinus++;
				}
			}
      }
      if (SIs[i].MatchedD == Plus) {
         LeftScore += SIs[i].score;
         LeftS++;
			if (SIs[i].Unique)
			   LeftUNum++;
      }
      else {
         RightScore += SIs[i].score;
         RightS++;
			if (SIs[i].Unique)
			   RightUNum++;
      }
   }
	
	
	short NumberSupSamples = 0; 
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
   for (short i = 0; i < VectorTag.size(); ++i) {
		if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) ++NumberSupSamples;
		if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) ++NumU_SupSamples;
		Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
	}
	//  << "\t" << Num_U_Reads
	//  << "\t" << NumU_SupSamples
	
   unsigned int EasyScore = LeftS * RightS;
   //double PreciseScore = (LeftScore + RightScore) * (-1);
   string CurrentReadSeq;
	
	CurrentChrMask[SIs[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[SIs[C_S].BPRight + SpacerBeforeAfter] = 'B';
	CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B';
	CurrentChrMask[RealEnd + SpacerBeforeAfter] = 'B';
	
   SIsOutf << "####################################################################################################" << endl;
   SIsOutf << NumberOfSIsInstances << "\tI " << SIs[C_S].IndelSize << "\tNT " << SIs[C_S].IndelSize << " \"" << SIs[C_S].InsertedStr << "\""
	<< "\tChrID " << SIs[C_S].FragName 
	<< "\tBP " << SIs[C_S].BPLeft + 1 << "\t" << SIs[C_S].BPRight + 1 
	<< "\tBP_range " << RealStart + 1 << "\t" << RealEnd + 1 
	<< "\tSupports " << NumberOfReads  << "\t" << Num_U_Reads
	<< "\t+ " << LeftS - 1  << "\t" << LeftUNum
	<< "\t- " << RightS - 1  << "\t" << RightUNum
	<< "\tS1 " << EasyScore; //EWL070111  << "\tS2 " << PreciseScore; 
	
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += SIs[i].MS;
   SIsOutf << "\tSUM_MS " << SUM_MS;
	
   SIsOutf << "\t" << VectorTag.size() << "\tNumSupSamples " << NumberSupSamples  << "\t" << NumU_SupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		SIsOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
		<< " " << NumSupportPerTag[i].NumUPlus
		<< " " << NumSupportPerTag[i].NumMinus
		<< " " << NumSupportPerTag[i].NumUMinus;
   SIsOutf << endl;
	
   SIsOutf << TheInput.substr(SIs[C_S].Left - ReportLength + SIs[C_S].BP + 1, ReportLength);// ReportLength
   for (int i = 0; i < SIs[C_S].IndelSize; i++) SIsOutf << " ";
   SIsOutf << TheInput.substr(SIs[C_S].Left + SIs[C_S].BP + 1, ReportLength) << endl;// ReportLength
   
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) {
      short SpaceBeforeReadSeq = ReportLength - SIs[GoodIndex].BP - 1;
      for (short i = 0; i < SpaceBeforeReadSeq; i++) SIsOutf << " "; 
      if (SIs[GoodIndex].MatchedD == Minus)
         SIsOutf << SIs[GoodIndex].UnmatchedSeq;
      else SIsOutf << ReverseComplement(SIs[GoodIndex].UnmatchedSeq);
      short SpaceBeforeD = ReportLength + ReportLength - SpaceBeforeReadSeq - SIs[GoodIndex].ReadLength;
      for (short i = 0; i < SpaceBeforeD; i++) SIsOutf << " "; 
      SIsOutf << "\t" << SIs[GoodIndex].MatchedD 
		<< "\t" << SIs[GoodIndex].MatchedRelPos 
		<< "\t" << SIs[GoodIndex].MS 
		<< "\t" << SIs[GoodIndex].Tag << "\t" << SIs[GoodIndex].Name << endl;  
   }
}

void OutputDI(const vector <SPLIT_READ> & DI,
				  const string & TheInput, 
				  const unsigned int & C_S,
				  const unsigned int & C_E,
				  const unsigned int & RealStart,
				  const unsigned int & RealEnd,
				  ofstream & DeletionOutf) {
   //short ReadLength = DI[C_S].ReadLength;
   //short ReadLengthMinus = ReadLength - 1;
   unsigned int NumberOfReads = C_E - C_S + 1;
   float LeftScore = 0;
   float RightScore = 0;
   unsigned int LeftS = 1;
   unsigned int RightS = 1;
	//int LeftNum = 0;
	int LeftUNum = 0;
	//int RightNum = 0;
	int RightUNum = 0;
	
   SupportPerSample NumSupportPerTag[VectorTag.size()];
	
   for (short i = 0; i < VectorTag.size(); i++) {
		NumSupportPerTag[i].NumPlus = 0;
		NumSupportPerTag[i].NumMinus = 0;
		NumSupportPerTag[i].NumUPlus = 0;
		NumSupportPerTag[i].NumUMinus = 0;
	}
	
   for (unsigned int i = C_S; i <= C_E; i++) {
      for (short j = 0; j < VectorTag.size(); j++) {
	      if (DI[i].Tag == VectorTag[j]) {
				if (DI[i].MatchedD == Plus) {
					NumSupportPerTag[j].NumPlus++;
					if (DI[i].Unique)
					   NumSupportPerTag[j].NumUPlus++;
				}
				else {
					NumSupportPerTag[j].NumMinus++;
					if (DI[i].Unique)
					   NumSupportPerTag[j].NumUMinus++;
				}
			}
      }
      if (DI[i].MatchedD == Plus) {
         LeftScore += DI[i].score;
         LeftS++;
			if (DI[i].Unique)
			   LeftUNum++;
      }
      else {
         RightScore += DI[i].score;
         RightS++;
			if (DI[i].Unique)
			   RightUNum++;
      }
   }
	
	short NumberSupSamples = 0; 
	short NumU_SupSamples = 0;
	int Num_U_Reads = 0;
   for (short i = 0; i < VectorTag.size(); ++i) {
		if (NumSupportPerTag[i].NumPlus + NumSupportPerTag[i].NumMinus) ++NumberSupSamples;
		if (NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus) ++NumU_SupSamples;
		Num_U_Reads += NumSupportPerTag[i].NumUPlus + NumSupportPerTag[i].NumUMinus;
	}
	//  << "\t" << Num_U_Reads
	
   unsigned int EasyScore = LeftS * RightS;
   //double PreciseScore = (LeftScore + RightScore) * (-1);
   //short GapSize =0;
   //if (DI[C_S].IndelSize < 14) GapSize = DI[C_S].IndelSize;
   //else GapSize = 13 + (int)log10(Deletions[C_S].IndelSize - 10);
	CurrentChrMask[DI[C_S].BPLeft + SpacerBeforeAfter] = 'B';
	CurrentChrMask[DI[C_S].BPRight + SpacerBeforeAfter] = 'B';
	//CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B';
	//CurrentChrMask[RealStart + SpacerBeforeAfter] = 'B'; 
   DeletionOutf << "####################################################################################################" << endl;
   DeletionOutf << NumberOfDIInstances + NumberOfDeletionsInstances + 1 << "\tD " << DI[C_S].IndelSize 
	<< "\tNT " << DI[C_S].NT_size << " \"" << DI[C_S].NT_str << "\"" 
	<< "\tChrID " << DI[C_S].FragName 
	<< "\tBP " << DI[C_S].BPLeft + 1 << "\t" << DI[C_S].BPRight + 1 
	<< "\tBP_range " << DI[C_S].BPLeft + 1 << "\t" << DI[C_S].BPRight + 1 
	<< "\tSupports " << NumberOfReads  << "\t" << Num_U_Reads
	<< "\t+ " << LeftS - 1  << "\t" << LeftUNum
	<< "\t- " << RightS - 1  << "\t" << RightUNum  
	<< "\tS1 " << EasyScore; // << "\tS2 0.0";// << PreciseScore << "\t"; 
	//EWL070111 << "\tS2 " << PreciseScore; 
	
   int SUM_MS = 0;
   for (unsigned int i = C_S; i <= C_E; i++) SUM_MS += DI[i].MS;
   DeletionOutf << "\tSUM_MS " << SUM_MS;
	
	
   DeletionOutf << "\t"<< VectorTag.size() << "\tNumSupSamples " << NumberSupSamples << "\t" << NumU_SupSamples;
   for (short i = 0; i < VectorTag.size(); i++)
		DeletionOutf << "\t" << VectorTag[i] << " " << NumSupportPerTag[i].NumPlus
		<< " " << NumSupportPerTag[i].NumUPlus
		<< " " << NumSupportPerTag[i].NumMinus
		<< " " << NumSupportPerTag[i].NumUMinus;
   DeletionOutf << endl;
	
   //DeletionOutf << TheInput.substr(DI[C_S].Left - ReportLength + DI[C_S].BP + 1, 2 * ReportLength) << endl;       
   DeletionOutf << TheInput.substr(DI[C_S].Left - ReportLength + DI[C_S].BP + 1, ReportLength);// << endl;// ReportLength    
	
   for (short i = 0; i < DI[C_S].NT_size; i++) DeletionOutf << " ";           
   DeletionOutf << TheInput.substr(DI[C_S].Left + DI[C_S].BP + 1 + DI[C_S].IndelSize, ReportLength) << endl;// ReportLength
   short SpaceBeforeReadSeq;    
   for (unsigned int GoodIndex = C_S; GoodIndex <= C_E; GoodIndex++) { 
      SpaceBeforeReadSeq = ReportLength - DI[GoodIndex].BP - 1;
      
      for (int i = 0; i < SpaceBeforeReadSeq; i++) DeletionOutf << " ";
      //short SpaceBeforeD = ReportLength + ReportLength - SpaceBeforeReadSeq - Deletions[GoodIndex].ReadLength;
      if (DI[GoodIndex].MatchedD == Minus) {
         DeletionOutf << DI[GoodIndex].UnmatchedSeq << "\t";
         //for (int i = 0; i < GapSize; i++) DeletionOutf << " ";    
         //DeletionOutf << Deletions[GoodIndex].UnmatchedSeq.substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      }
      else {
         DeletionOutf << ReverseComplement(DI[GoodIndex].UnmatchedSeq) << "\t"; 
         //for (int i = 0; i < GapSize; i++) DeletionOutf << " ";  
         //DeletionOutf << ReverseComplement(Deletions[GoodIndex].UnmatchedSeq).substr(Deletions[GoodIndex].BP + 1, Deletions[GoodIndex].ReadLength - Deletions[GoodIndex].BP);// << endl;
      } 
      //for (int i = 0; i < SpaceBeforeD; i++) DeletionOutf << " ";
      DeletionOutf << "\t" << DI[GoodIndex].MatchedD << "\t" << DI[GoodIndex].MatchedRelPos
		<< "\t" << DI[GoodIndex].MS
		<< "\t" << DI[GoodIndex].Tag << "\t" <<  DI[GoodIndex].Name << endl;
   }     
}

void SortOutputSI(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> SIs[], ofstream & SIsOutf) {
   cout << "Sorting and outputing short insertions ..." << endl;
   unsigned int SIsNum;
   short CompareResult;
   unsigned Temp4Exchange;
	
	//vector <SPLIT_READ> InputIndels; 
   vector <SPLIT_READ> GoodIndels;
   unsigned int GoodNum;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (SIs[Box_index].size() >= NumRead2ReportCutOff) {
			//InputIndels.clear();
			SIsNum = SIs[Box_index].size();
			//cout << "SIsNum " << SIsNum << endl;
			//for (int i = 0; i < SIsNum; i++) {
			//	InputIndels.push_back(Reads[SIs[Box_index][i]]);
			//}
			for (unsigned int First = 0; First < SIsNum - 1; First++) {
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < SIsNum; Second++) {
                  //if (InputIndels[Second].Unique) 
						{  
							if (Reads[SIs[Box_index][First]].ReadLength == Reads[SIs[Box_index][Second]].ReadLength) {
								if (Reads[SIs[Box_index][First]].LeftMostPos == Reads[SIs[Box_index][Second]].LeftMostPos)
									Reads[SIs[Box_index][Second]].Unique = false;
							}
							
							if (Reads[SIs[Box_index][First]].BPLeft < Reads[SIs[Box_index][Second]].BPLeft) continue;
							else if (Reads[SIs[Box_index][First]].BPLeft > Reads[SIs[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else {
								if (Reads[SIs[Box_index][First]].IndelSize < Reads[SIs[Box_index][Second]].IndelSize) continue;
								else if (Reads[SIs[Box_index][First]].IndelSize > Reads[SIs[Box_index][Second]].IndelSize) {
									CompareResult = 1;
								}
								//else { // InputIndels[First].BPRight == InputIndels[Second].BPRight
								//		short Compare2Str = CompareTwoString(InputIndels[First].InsertedStr, InputIndels[Second].InsertedStr );
								//if (Compare2Str > 0) CompareResult = 1;
								//else if (Compare2Str == 0) CompareResult = 2;
								//else continue;
								
								//}
								//else if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//	if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//		InputIndels[Second].Unique = false;
								//	}
								//}
							}
                     //CompareResult = CompareTwoReads(InputIndels[First], InputIndels[Second]);
                     if (CompareResult == 1) {
  	                     Temp4Exchange = SIs[Box_index][First];
  	                     SIs[Box_index][First] = SIs[Box_index][Second];
  	                     SIs[Box_index][Second] = Temp4Exchange;
                     }
                     //else if (CompareResult == 2) InputIndels[Second].Unique = false;
                  }
               }
            }
         }
         GoodIndels.clear();
	      IndelEvents.clear();
			//cout << GoodIndels.size() << endl;
			
			for (unsigned int First = 0; First < SIsNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(Reads[SIs[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size();
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
         Indel4output OneIndelEvent;
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.IndelSize = GoodIndels[0].IndelSize;
         OneIndelEvent.IndelStr = GoodIndels[0].InsertedStr;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
			OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
         //OneIndelEvent.IndelStr = GoodIndels[0].InsertedStr;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].IndelSize == OneIndelEvent.IndelSize)
					//&& OneIndelEvent.IndelStr == GoodIndels[GoodIndex].InsertedStr )
					
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					GetRealStart4Insertion(CurrentChr, OneIndelEvent.IndelStr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
					OneIndelEvent.IndelSize =  GoodIndels[GoodIndex].IndelSize;
					OneIndelEvent.IndelStr =  GoodIndels[GoodIndex].InsertedStr;
					//OneIndelEvent.IndelStr = GoodIndels[GoodIndex].InsertedStr;
            }
         }
			
			//if (OneIndelEvent.End - OneIndelEvent.Start + 1 >= NumRead2ReportCutOff) 
			OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
			OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
			
			OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			GetRealStart4Insertion(CurrentChr, OneIndelEvent.IndelStr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
			IndelEvents.push_back(OneIndelEvent);
			
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			unsigned int IndelSize;
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  IndelSize = IndelEvents[EventIndex].IndelSize;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;
						
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else if (IndelEvents[EventIndex_left].IndelSize != IndelSize) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        } 
                     }
                  }
                  // report max one
                  //cout << Max_Support << endl;
                  if (Max_Support >= NumRead2ReportCutOff) {
							OutputSIs(GoodIndels, CurrentChr,
										 IndelEvents[Max_Support_Index].Start,
										 IndelEvents[Max_Support_Index].End,
										 RealStart, RealEnd, SIsOutf);
							NumberOfSIsInstances++;
						}
               }
            }
         }
      }   // if (!insertion[Box_index].empty())
   }
   cout << "Short insertions: " << NumberOfSIsInstances << endl << endl;
}

void SortOutputTD(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> TDs[], ofstream & TDOutf) {
   cout << "Sorting and outputing tandem duplications ..." << endl;
   unsigned int TDNum;
   short CompareResult;
   unsigned Temp4Exchange;
	
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (TDs[Box_index].size() >= NumRead2ReportCutOff) {
			//InputIndels.clear();
			TDNum = TDs[Box_index].size();
			//for (int i = 0 ; i < TDNum; i ++) {
			//	InputIndels.push_back(AllReads[TDs[Box_index][i]]);
			//}
         
         for (unsigned int First = 0; First < TDNum - 1; First++) {
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < TDNum; Second++) {
                  //if (InputIndels[Second].Unique) 
						{
							if (AllReads[TDs[Box_index][First]].ReadLength == AllReads[TDs[Box_index][Second]].ReadLength) {
								if (AllReads[TDs[Box_index][First]].LeftMostPos == AllReads[TDs[Box_index][Second]].LeftMostPos)
									AllReads[TDs[Box_index][Second]].Unique = false;
							}
							if (AllReads[TDs[Box_index][First]].BPLeft < AllReads[TDs[Box_index][Second]].BPLeft) continue;
                     else if (AllReads[TDs[Box_index][First]].BPLeft > AllReads[TDs[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (AllReads[TDs[Box_index][First]].BPLeft == AllReads[TDs[Box_index][Second]].BPLeft) {
								if (AllReads[TDs[Box_index][First]].BPRight < AllReads[TDs[Box_index][Second]].BPRight) continue;
								else if (AllReads[TDs[Box_index][First]].BPRight > AllReads[TDs[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].Unique = false;
								//		}
								//		
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = TDs[Box_index][First];
  	                     TDs[Box_index][First] = TDs[Box_index][Second];
  	                     TDs[Box_index][Second] = Temp4Exchange;
							}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();
			
			for (unsigned int First = 0; First < TDNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(AllReads[TDs[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size(); 
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft 
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }
			
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//	     cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        } 
                     }
                  }
                  // report max one
                  if (IndelEvents[Max_Support_Index].Support >= NumRead2ReportCutOff) {
							if (GoodIndels[IndelEvents[Max_Support_Index].Start].IndelSize < BalanceCutoff) {
								OutputTDs(GoodIndels, CurrentChr, 
											 IndelEvents[Max_Support_Index].Start, 
											 IndelEvents[Max_Support_Index].End, 
											 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
							}
							else if (ReportEvent(GoodIndels, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End)) {
								OutputTDs(GoodIndels, CurrentChr, 
											 IndelEvents[Max_Support_Index].Start, 
											 IndelEvents[Max_Support_Index].End, 
											 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
							}
						}
               }
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Tandem duplications: " << NumberOfTDInstances << endl << endl;
}

void SortOutputTD_NT(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & AllReads, vector <unsigned> TDs[], ofstream & TDOutf) {
   cout << "Sorting and outputing tandem duplications with non-template sequence ..." << endl;
   unsigned int TDNum;
   short CompareResult;
   unsigned Temp4Exchange;
	
	int Count_TD_NT_output = 0;
	
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      if (TDs[Box_index].size() >= NumRead2ReportCutOff) {
			//InputIndels.clear();
			TDNum = TDs[Box_index].size();
			//for (int i = 0 ; i < TDNum; i ++) {
			//	InputIndels.push_back(AllReads[TDs[Box_index][i]]);
			//}
         
         for (unsigned int First = 0; First < TDNum - 1; First++) {
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < TDNum; Second++) {
                  //if (InputIndels[Second].Unique) 
						{
							if (AllReads[TDs[Box_index][First]].ReadLength == AllReads[TDs[Box_index][Second]].ReadLength) {
								if (AllReads[TDs[Box_index][First]].LeftMostPos == AllReads[TDs[Box_index][Second]].LeftMostPos)
									AllReads[TDs[Box_index][Second]].Unique = false;
							}
							if (AllReads[TDs[Box_index][First]].BPLeft < AllReads[TDs[Box_index][Second]].BPLeft) continue;
                     else if (AllReads[TDs[Box_index][First]].BPLeft > AllReads[TDs[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (AllReads[TDs[Box_index][First]].BPLeft == AllReads[TDs[Box_index][Second]].BPLeft) {
								if (AllReads[TDs[Box_index][First]].BPRight < AllReads[TDs[Box_index][Second]].BPRight) continue;
								else if (AllReads[TDs[Box_index][First]].BPRight > AllReads[TDs[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								else { // InputIndels[First].BPRight == InputIndels[Second].BPRight
									if (AllReads[TDs[Box_index][First]].NT_size < AllReads[TDs[Box_index][Second]].NT_size) continue;
									else if (AllReads[TDs[Box_index][First]].NT_size > AllReads[TDs[Box_index][Second]].NT_size) CompareResult = 1;
									//else { // InputIndels[First].NT_size == InputIndels[Second].NT_size
									//	short Compare2Str = CompareTwoString(InputIndels[First].NT_str, InputIndels[Second].NT_str );
									//	if (Compare2Str > 0) CompareResult = 1;
									//	else if (Compare2Str == 0) CompareResult = 2;
									//	else continue;
									//}
								}
							}
							if (CompareResult == 1) {
								Temp4Exchange = TDs[Box_index][First];
  	                     TDs[Box_index][First] = TDs[Box_index][Second];
  	                     TDs[Box_index][Second] = Temp4Exchange;
							}
							//else if (CompareResult == 2) {
							//	Temp4Exchange = InputIndels[First + 1];
							//  InputIndels[First + 1] = InputIndels[Second];
							// InputIndels[Second] = Temp4Exchange;
							//}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();
			
			for (unsigned int First = 0; First < TDNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(AllReads[TDs[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size(); 
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft 
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }
			
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//	     cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        } 
                     }
                  }
                  // report max one
                  if (IndelEvents[Max_Support_Index].Support >= NumRead2ReportCutOff) {
							if (GoodIndels[IndelEvents[Max_Support_Index].Start].IndelSize < BalanceCutoff) {
								OutputTDs(GoodIndels, CurrentChr, 
											 IndelEvents[Max_Support_Index].Start, 
											 IndelEvents[Max_Support_Index].End, 
											 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
								Count_TD_NT_output++;
							}
							else if (ReportEvent(GoodIndels, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End)) {
								OutputTDs(GoodIndels, CurrentChr, 
											 IndelEvents[Max_Support_Index].Start, 
											 IndelEvents[Max_Support_Index].End, 
											 RealStart, RealEnd, TDOutf);
								NumberOfTDInstances++;
								Count_TD_NT_output++;
							}
						}
               }
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Tandem duplications with non-template sequence (TD_NT): " << Count_TD_NT_output << "\n\n";
}

void SortOutputD(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> Deletions[], ofstream & DeletionOutf) {
   cout << "Sorting and outputing deletions ..." << endl;
   unsigned int DeletionsNum;
   short CompareResult;
   unsigned Temp4Exchange;
	
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
		//cout << Box_index << "\t" << NumBoxes << "\t" << Deletions[Box_index].size() << endl;
		//if (Deletions[Box_index].size() >= NumRead2ReportCutOff)
		//   cout << Box_index << "\t" << Deletions[Box_index].size() << endl;
      if (Deletions[Box_index].size() >= NumRead2ReportCutOff) {
			DeletionsNum = Deletions[Box_index].size();
         
         for (unsigned int First = 0; First < DeletionsNum - 1; First++) {
            //if (InputIndels[First].Unique) 
				//cout << DeletionsNum << " First " << First << endl;
				{
               for (unsigned int Second = First + 1; Second < DeletionsNum; Second++) {
                  //if (InputIndels[Second].Unique) 
						//cout << DeletionsNum << " Second " << Second << endl;
						{
							if (Reads[Deletions[Box_index][First]].ReadLength == Reads[Deletions[Box_index][Second]].ReadLength) {
								if (Reads[Deletions[Box_index][First]].LeftMostPos == Reads[Deletions[Box_index][Second]].LeftMostPos)
									Reads[Deletions[Box_index][Second]].Unique = false;
							}
							if (Reads[Deletions[Box_index][First]].BPLeft < Reads[Deletions[Box_index][Second]].BPLeft) continue;
                     else if (Reads[Deletions[Box_index][First]].BPLeft > Reads[Deletions[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (Reads[Deletions[Box_index][First]].BPLeft == Reads[Deletions[Box_index][Second]].BPLeft) {
								if (Reads[Deletions[Box_index][First]].BPRight < Reads[Deletions[Box_index][Second]].BPRight) continue;
								else if (Reads[Deletions[Box_index][First]].BPRight > Reads[Deletions[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								//else CompareResult = 2;
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].Unique = false;
								//		}
								//			
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = Deletions[Box_index][First];
  	                     Deletions[Box_index][First] = Deletions[Box_index][Second];
  	                     Deletions[Box_index][Second] = Temp4Exchange;
							}
							//else if (CompareResult == 2) {
							//	Temp4Exchange = InputIndels[First + 1];
  	                  //   InputIndels[First + 1] = InputIndels[Second];
  	                  //   InputIndels[Second] = Temp4Exchange;
							//}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();
			//for (int i = 0 ; i < DeletionsNum; i ++) {
			//	InputIndels.push_back(Reads[Deletions[Box_index][i]]);
			//}Reads[Deletions[Box_index][First]]
         //cout << "GoodIndels" << endl;
			for (unsigned int First = 0; First < DeletionsNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(Reads[Deletions[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size();
			//cout << Box_index << " box read size " << GoodNum << endl;
         if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
			//cout << "here" << endl;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft 
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }
			
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
				     //cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
					//cout << EventIndex << " EventIndex" << endl;
               if (IndelEvents[EventIndex].WhetherReport) {
                  RealStart = IndelEvents[EventIndex].RealStart;
                  RealEnd = IndelEvents[EventIndex].RealEnd;
                  Max_Support = IndelEvents[EventIndex].Support;
                  Max_Support_Index = EventIndex;
                  for (unsigned EventIndex_left = 0; EventIndex_left < IndelEvents.size(); EventIndex_left++) {
                     if (IndelEvents[EventIndex_left].WhetherReport == false) continue;
                     else if (IndelEvents[EventIndex_left].RealStart != RealStart) continue;
                     else if (IndelEvents[EventIndex_left].RealEnd != RealEnd) continue;
                     else {
                        IndelEvents[EventIndex_left].WhetherReport = false;
                        if (IndelEvents[EventIndex_left].Support > Max_Support) {
                           Max_Support = IndelEvents[EventIndex_left].Support;
                           Max_Support_Index = EventIndex_left;
                        } 
                     }
                  }
                  // report max one
						//cout << "max" << endl;
                  if (IndelEvents[Max_Support_Index].Support >= NumRead2ReportCutOff) {
							//cout << "aa" << endl;
							if (GoodIndels[IndelEvents[Max_Support_Index].Start].IndelSize < BalanceCutoff) {
								//cout << "ba" << endl;
								OutputDeletions(GoodIndels, CurrentChr, 
													 IndelEvents[Max_Support_Index].Start, 
													 IndelEvents[Max_Support_Index].End, 
													 RealStart, RealEnd, DeletionOutf);
								NumberOfDeletionsInstances++;
								//cout << "bb" << endl;
							}
							else if (ReportEvent(GoodIndels, IndelEvents[Max_Support_Index].Start, IndelEvents[Max_Support_Index].End)) {
								//cout << "ca" << endl;
								OutputDeletions(GoodIndels, CurrentChr, 
													 IndelEvents[Max_Support_Index].Start, 
													 IndelEvents[Max_Support_Index].End, 
													 RealStart, RealEnd, DeletionOutf);
								NumberOfDeletionsInstances++;
								//cout << "cb" << endl;
							}
						}
               }
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Deletions: " << NumberOfDeletionsInstances << endl << endl;
}

void SortOutputInv(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> Inv[], ofstream & InvOutf) {
   cout << "Sorting and outputing Inversions ..." << endl;
   unsigned int InversionsNum;
   short CompareResult;
   unsigned Temp4Exchange;
	/*
	 unsigned int C_S = 0;
	 unsigned int C_E = 0;
	 unsigned int C_BP_Left;// = GoodSIs[0].BPLeft;
	 unsigned int C_BP_Right;// = GoodSIs[0].BPRight;
	 unsigned int C_Indelsize;
	 */
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
		
      if (Inv[Box_index].size() >= NumRead2ReportCutOff) {
			//cout << Box_index << "\t" << Inv[Box_index].size() << endl;
         InversionsNum = Inv[Box_index].size();
			//InputIndels.clear();
			//for (int i = 0; i < InversionsNum; i++) InputIndels.push_back(Reads[Inv[Box_index][i]]);
         for (unsigned int First = 0; First < InversionsNum - 1; First++) {
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < InversionsNum; Second++) {
                  //if (InputIndels[Second].Unique) 
						{
							if (Reads[Inv[Box_index][First]].ReadLength == Reads[Inv[Box_index][Second]].ReadLength) {
								if (Reads[Inv[Box_index][First]].LeftMostPos == Reads[Inv[Box_index][Second]].LeftMostPos)
									Reads[Inv[Box_index][Second]].Unique = false;
							}
							if (Reads[Inv[Box_index][First]].BPLeft < Reads[Inv[Box_index][Second]].BPLeft) continue;
                     else if (Reads[Inv[Box_index][First]].BPLeft > Reads[Inv[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (Reads[Inv[Box_index][First]].BPLeft == Reads[Inv[Box_index][Second]].BPLeft) {
								if (Reads[Inv[Box_index][First]].BPRight < Reads[Inv[Box_index][Second]].BPRight) continue;
								else if (Reads[Inv[Box_index][First]].BPRight > Reads[Inv[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								//else CompareResult = 2;
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].Unique = false;
								//		}
								//			
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = Inv[Box_index][First];
  	                     Inv[Box_index][First] = Inv[Box_index][Second];
  	                     Inv[Box_index][Second] = Temp4Exchange;
							}
							//else if (CompareResult == 2) {
							//	Temp4Exchange = InputIndels[First + 1];
  	                  //   InputIndels[First + 1] = InputIndels[Second];
  	                  //   InputIndels[Second] = Temp4Exchange;
							//}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();
			
			for (unsigned int First = 0; First < InversionsNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(Reads[Inv[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size();  
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft 
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					//GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }
			
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
         //GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//	     cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
					RealStart = IndelEvents[EventIndex].RealStart;
					RealEnd = IndelEvents[EventIndex].RealEnd;
					if (IndelEvents[EventIndex].Support < NumRead2ReportCutOff) continue;
					// report max one
					if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < BalanceCutoff) {
						OutputInversions(GoodIndels, CurrentChr, 
											  IndelEvents[EventIndex].Start, 
											  IndelEvents[EventIndex].End, 
											  RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;						
					}
					else if (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) {
						OutputInversions(GoodIndels, CurrentChr, 
											  IndelEvents[EventIndex].Start, 
											  IndelEvents[EventIndex].End, 
											  RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;
					}
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Inversions (INV): " << NumberOfInvInstances << endl << endl;
}

void SortOutputInv_NT(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> Inv[], ofstream & InvOutf) {
   cout << "Sorting and outputing Inversions with non-template sequence ..." << endl;
   unsigned int InversionsNum;
   short CompareResult;
   unsigned Temp4Exchange;
	
	int Count_INV_NT_output = 0;
	/*
	 unsigned int C_S = 0;
	 unsigned int C_E = 0;
	 unsigned int C_BP_Left;// = GoodSIs[0].BPLeft;
	 unsigned int C_BP_Right;// = GoodSIs[0].BPRight;
	 unsigned int C_Indelsize;
	 */
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
		
      if (Inv[Box_index].size() >= NumRead2ReportCutOff) {
         InversionsNum = Inv[Box_index].size();
			//cout << Box_index << "\t" << Inv[Box_index].size() << endl;
			//InputIndels.clear();
			//for (int i = 0; i < InversionsNum; i++) InputIndels.push_back(Reads[Inv[Box_index][i]]);
         for (unsigned int First = 0; First < InversionsNum - 1; First++) {
				
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < InversionsNum; Second++) {
						//cout << InputIndels[First].BPLeft << "\t" << InputIndels[First].BPRight << "\t" 
						//     << InputIndels[Second].BPLeft << "\t" << InputIndels[Second].BPRight<< endl;
                  //if (InputIndels[Second].Unique) 
						{
							if (Reads[Inv[Box_index][First]].ReadLength == Reads[Inv[Box_index][Second]].ReadLength) {
								if (Reads[Inv[Box_index][First]].LeftMostPos == Reads[Inv[Box_index][Second]].LeftMostPos)
									Reads[Inv[Box_index][Second]].Unique = false;
							}
							if (Reads[Inv[Box_index][First]].BPLeft < Reads[Inv[Box_index][Second]].BPLeft) continue;
                     else if (Reads[Inv[Box_index][First]].BPLeft > Reads[Inv[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (Reads[Inv[Box_index][First]].BPLeft == Reads[Inv[Box_index][Second]].BPLeft) {
								if (Reads[Inv[Box_index][First]].BPRight < Reads[Inv[Box_index][Second]].BPRight) continue;
								else if (Reads[Inv[Box_index][First]].BPRight > Reads[Inv[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								else { // InputIndels[First].BPRight == InputIndels[Second].BPRight
									if (Reads[Inv[Box_index][First]].NT_size < Reads[Inv[Box_index][Second]].NT_size) continue;
									else if (Reads[Inv[Box_index][First]].NT_size > Reads[Inv[Box_index][Second]].NT_size) CompareResult = 1;
									//else { // InputIndels[First].NT_size == InputIndels[Second].NT_size
									//	short Compare2Str = CompareTwoString(InputIndels[First].NT_str, InputIndels[Second].NT_str );
									//	if (Compare2Str > 0) CompareResult = 1;
									//	else if (Compare2Str == 0) CompareResult = 2;
									//	else continue;
									//}
								}
								//else {
								//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
								//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
								//			InputIndels[Second].Unique = false;
								//		}
								//		
								//	}
								//}
							}
							if (CompareResult == 1) {
								Temp4Exchange = Inv[Box_index][First];
  	                     Inv[Box_index][First] = Inv[Box_index][Second];
  	                     Inv[Box_index][Second] = Temp4Exchange;
							}
							//if (CompareResult == 2) {
							//	Temp4Exchange = InputIndels[First + 1];
  	                  //   InputIndels[First + 1] = InputIndels[Second];
  	                  //   InputIndels[Second] = Temp4Exchange;
							//}
                  }
               }
            }
         }
         GoodIndels.clear();
			IndelEvents.clear();
			
			for (unsigned int First = 0; First < InversionsNum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(Reads[Inv[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size();  
			//cout << "GoodNum " << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
			//    cout << GoodNum << endl;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
         OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
				//cout << GoodIndex << "\t" << GoodIndels[GoodIndex].BPLeft << "\t" << GoodIndels[GoodIndex].BPRight << "\t" << OneIndelEvent.BPLeft << "\t" << OneIndelEvent.BPRight << endl;
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft 
                && GoodIndels[GoodIndex].BPRight == OneIndelEvent.BPRight)
               OneIndelEvent.End = GoodIndex;
				else  {
					OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
					OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
					OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
					//GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.BPRight =  GoodIndels[GoodIndex].BPRight;
            }
         }
			
         OneIndelEvent.RealStart = OneIndelEvent.BPLeft;
         OneIndelEvent.RealEnd = OneIndelEvent.BPRight;
         OneIndelEvent.Support = OneIndelEvent.End - OneIndelEvent.Start + 1;
			//cout << OneIndelEvent.Support << "\t" << OneIndelEvent.Start << "\t" << OneIndelEvent.End << endl;
         //GetRealStart4Deletion(CurrentChr, OneIndelEvent.RealStart, OneIndelEvent.RealEnd);
         IndelEvents.push_back(OneIndelEvent);
			//cout << "IndelEvent: " << IndelEvents.size() << endl;
			//string IndelType;
			unsigned int RealStart;
			unsigned int RealEnd;
			//bool WhetherDeletion = true;
			string IndelStr;
			unsigned int Max_Support;
			unsigned int Max_Support_Index;
			
			if (IndelEvents.size()) {
            for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
					//cout << IndelEvents[EventIndex].Start << "\t" << IndelEvents[EventIndex].End << "\t" << IndelEvents[EventIndex].Support << endl;
					RealStart = IndelEvents[EventIndex].RealStart;
					RealEnd = IndelEvents[EventIndex].RealEnd;
					if (IndelEvents[EventIndex].Support < NumRead2ReportCutOff) continue;
					// report max one
					if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < BalanceCutoff) {
						OutputInversions(GoodIndels, CurrentChr, 
											  IndelEvents[EventIndex].Start, 
											  IndelEvents[EventIndex].End, 
											  RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;
						Count_INV_NT_output++;
					}
					else if (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) {
						OutputInversions(GoodIndels, CurrentChr, 
											  IndelEvents[EventIndex].Start, 
											  IndelEvents[EventIndex].End, 
											  RealStart, RealEnd, InvOutf);
						NumberOfInvInstances++;
						Count_INV_NT_output++;
					}
            }
         }
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "Inversions with non-template sequence (INV_NT): " << Count_INV_NT_output << endl << endl;
}

void SortOutputDI(const unsigned & NumBoxes, const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <unsigned> DI[], ofstream & DIOutf) {
   cout << "Sorting and outputing deletions with non-template sequences ..." << endl;
   unsigned int DINum;
   short CompareResult;
   unsigned Temp4Exchange;
	/*
	 unsigned int C_S = 0;
	 unsigned int C_E = 0;
	 unsigned int C_BP_Left;// = GoodSIs[0].BPLeft;
	 unsigned int C_BP_Right;// = GoodSIs[0].BPRight;
	 unsigned int C_Indelsize;
	 */
   unsigned int GoodNum;
	//vector <SPLIT_READ> InputIndels;
   vector <SPLIT_READ> GoodIndels;
   vector <Indel4output> IndelEvents;
	
   for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++) {
      //cout << "Box_index: "   << Box_index << endl;
      if (DI[Box_index].size() >= NumRead2ReportCutOff) {
         DINum = DI[Box_index].size();
			//InputIndels.clear();
			//for (int i = 0; i < DINum; i++) InputIndels.push_back(Reads[DI[Box_index][i]]);
         for (unsigned int First = 0; First < DINum - 1; First++) {
            //if (InputIndels[First].Unique) 
				{
               for (unsigned int Second = First + 1; Second < DINum; Second++) {
                  //if (InputIndels[Second].Unique) 
						{
							if (Reads[DI[Box_index][First]].ReadLength == Reads[DI[Box_index][Second]].ReadLength) {
								if (Reads[DI[Box_index][First]].LeftMostPos == Reads[DI[Box_index][Second]].LeftMostPos)
									Reads[DI[Box_index][Second]].Unique = false;
							}
							if (Reads[DI[Box_index][First]].BPLeft < Reads[DI[Box_index][Second]].BPLeft) continue;
                     else if (Reads[DI[Box_index][First]].BPLeft > Reads[DI[Box_index][Second]].BPLeft) {
								CompareResult = 1;
							}
							else if (Reads[DI[Box_index][First]].BPLeft == Reads[DI[Box_index][Second]].BPLeft) {
								if (Reads[DI[Box_index][First]].BPRight 
									 < Reads[DI[Box_index][Second]].BPRight) continue;
								else if (Reads[DI[Box_index][First]].BPRight > Reads[DI[Box_index][Second]].BPRight) {
									CompareResult = 1;
								}
								else {
									if (Reads[DI[Box_index][First]].NT_size < Reads[DI[Box_index][Second]].NT_size) continue;
									else if (Reads[DI[Box_index][First]].NT_size > Reads[DI[Box_index][Second]].NT_size) CompareResult = 1;
									//else CompareResult = 2;
									//else {
									//	if (InputIndels[First].MatchedRelPos == InputIndels[Second].MatchedRelPos) {
									//		if (InputIndels[First].UnmatchedSeq == InputIndels[Second].UnmatchedSeq) {
									//			InputIndels[Second].Unique = false;
									//		}
									//			
									//	}
									//}
							   }
                     }
							if (CompareResult == 1) {
								Temp4Exchange = DI[Box_index][First];
  	                     DI[Box_index][First] = DI[Box_index][Second];
  	                     DI[Box_index][Second] = Temp4Exchange;
							}
							//else if (CompareResult == 2) {
							//	Temp4Exchange = InputIndels[First + 1];
  	                  //   InputIndels[First + 1] = InputIndels[Second];
  	                  //   InputIndels[Second] = Temp4Exchange;
							//}
						}
               }
            }
         }
         GoodIndels.clear();
	      IndelEvents.clear();
			
			for (unsigned int First = 0; First < DINum; First++) {
				//if (InputIndels[First].Unique) 
				GoodIndels.push_back(Reads[DI[Box_index][First]]);
			}				
			
         GoodNum = GoodIndels.size(); 
			//cout << Box_index << " " << GoodNum << endl;
			if (GoodNum == 0) continue;
         //cout << "GoodNum: " << GoodNum << endl;   string InsertedStr;   string NT_str;  short NT_size;
         Indel4output OneIndelEvent; 
         OneIndelEvent.Start = 0;
			OneIndelEvent.End = 0;
			OneIndelEvent.IndelSize = GoodIndels[0].IndelSize;
         OneIndelEvent.NT_size = GoodIndels[0].NT_size;
			OneIndelEvent.BPLeft =  GoodIndels[0].BPLeft;
			OneIndelEvent.BPRight =  GoodIndels[0].BPRight;
			//OneIndelEvent.IndelStr = GoodIndels[0].NT_str;
			OneIndelEvent.WhetherReport = true;
         for (unsigned int GoodIndex = 1; GoodIndex < GoodNum; GoodIndex++) {
            if (GoodIndels[GoodIndex].BPLeft == OneIndelEvent.BPLeft
                && GoodIndels[GoodIndex].IndelSize == OneIndelEvent.IndelSize  
                && GoodIndels[GoodIndex].NT_size == OneIndelEvent.NT_size)
					//&& OneIndelEvent.IndelStr == GoodIndels[GoodIndex].NT_str)
               OneIndelEvent.End = GoodIndex;
				else  {
					IndelEvents.push_back(OneIndelEvent);
					OneIndelEvent.Start = GoodIndex;
					OneIndelEvent.End = GoodIndex;
					OneIndelEvent.BPLeft =  GoodIndels[GoodIndex].BPLeft;
					OneIndelEvent.IndelSize =  GoodIndels[GoodIndex].IndelSize;
					OneIndelEvent.NT_size = GoodIndels[GoodIndex].NT_size;
					//OneIndelEvent.IndelStr = GoodIndels[GoodIndex].NT_str;
            }
         }
			
			//if (OneIndelEvent.End - OneIndelEvent.Start + 1 >= NumRead2ReportCutOff) 
			IndelEvents.push_back(OneIndelEvent);
			unsigned int RealStart;
			unsigned int RealEnd;
			for (unsigned EventIndex = 0; EventIndex < IndelEvents.size(); EventIndex++) {
				if (IndelEvents[EventIndex].End - IndelEvents[EventIndex].Start + 1 >= NumRead2ReportCutOff) 
				{
					RealStart = GoodIndels[IndelEvents[EventIndex].Start].BPLeft;
					RealEnd = GoodIndels[IndelEvents[EventIndex].Start].BPRight;
					//if (IndelEvents[EventIndex].IndelSize < 100) {}
					//if (IndelSize < BalanceCutoff) {
					//        OutputDI(GoodIndels, CurrentChr, 
					//                 IndelEvents[EventIndex].Start, 
					//       	        IndelEvents[EventIndex].End, 
					// 		RealStart, RealEnd, DIOutf);
					//        NumberOfDIInstances++;
					//}
					if (GoodIndels[IndelEvents[EventIndex].Start].IndelSize < BalanceCutoff) {
						OutputDI(GoodIndels, CurrentChr, 
									IndelEvents[EventIndex].Start, 
									IndelEvents[EventIndex].End, 
									RealStart, RealEnd, DIOutf);
						NumberOfDIInstances++;					
					}
					else if (ReportEvent(GoodIndels, IndelEvents[EventIndex].Start, IndelEvents[EventIndex].End)) {
						OutputDI(GoodIndels, CurrentChr, 
									IndelEvents[EventIndex].Start, 
									IndelEvents[EventIndex].End, 
									RealStart, RealEnd, DIOutf);
						NumberOfDIInstances++;
					}
				}
			} 
      }   // if (!Deletions[Box_index].empty())
   } // for (unsigned Box_index = 0; Box_index < NumBoxes; Box_index++)
   cout << "deletions with non-template sequences: " << NumberOfDIInstances << endl << endl;
}

void SortOutputLI(const string & CurrentChr, vector <SPLIT_READ> & Reads, ofstream & LargeInsertionOutf) {
	unsigned UP_Close_index;
	unsigned temp_AbsLoc;
	// find LI combinations
	uint8_t * plus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	uint8_t * minus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	int32_t * EventIndex_Pos = new int32_t[CurrentChr.size() + 1];
	for (unsigned i = 0; i < CurrentChr.size() + 1; i++) {
		plus_LI_Pos[i] = 0;
		minus_LI_Pos[i] = 0;
		EventIndex_Pos[i] = -1;
	}
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
		//UP_Close_index = ;
		if (Reads[Index].Found || Reads[Index].Used || !Reads[Index].UP_Far.empty()) continue;
		temp_AbsLoc = Reads[Index].UP_Close[Reads[Index].UP_Close.size() - 1].AbsLoc;
		if (plus_LI_Pos[temp_AbsLoc] < Max_short) {
			if (Reads[Index].MatchedD == Plus) plus_LI_Pos[temp_AbsLoc]++;	
		}
		if (minus_LI_Pos[temp_AbsLoc] < Max_short) {
			if (Reads[Index].MatchedD == Minus) minus_LI_Pos[temp_AbsLoc]++;
		}
	}
	vector <LI_Pos> LI_Positions;
	LI_Pos temp_LI_pos;
	bool SkipThis;
	int LI_Positions_Size = 0;
	bool SkipPlus;
	for (int Index_Minus = SpacerBeforeAfter; Index_Minus < CurrentChr.size() - SpacerBeforeAfter; Index_Minus++) {
		SkipPlus = false;
		for (int MaskedPosIndexMinus = Index_Minus + 10; MaskedPosIndexMinus >= Index_Minus - 10; MaskedPosIndexMinus--) {
			if (CurrentChrMask[MaskedPosIndexMinus] == 'B') {
				Index_Minus = MaskedPosIndexMinus + 10;
				SkipPlus = true;
				break;
			} 
		}
		if (SkipPlus == false && minus_LI_Pos[Index_Minus] >= NumRead2ReportCutOff) {
			for (int Index_Plus = Index_Minus - 1; Index_Plus <= Index_Minus + 30; Index_Plus++) {
				SkipThis = false;
				for (int MaskedPosIndexPlus = Index_Plus + 10; MaskedPosIndexPlus >= Index_Plus - 10; MaskedPosIndexPlus--) {
					if (CurrentChrMask[MaskedPosIndexPlus] == 'B') {
						if ((MaskedPosIndexPlus + 10) > Index_Minus ) { Index_Minus = MaskedPosIndexPlus + 10; }
						SkipThis = true;
						break;
					} 
				}
				if (SkipThis == false && plus_LI_Pos[Index_Plus] >= NumRead2ReportCutOff) {
					temp_LI_pos.Plus_Pos = Index_Plus;
					temp_LI_pos.Minus_Pos = Index_Minus;
					CurrentChrMask[Index_Plus] == 'B';
					CurrentChrMask[Index_Minus] == 'B';
					//cout << Index_Plus << "\t" << (short)plus_LI_Pos[Index_Plus] << "\t" 
					//     << Index_Minus << "\t" << (short)minus_LI_Pos[Index_Minus] << endl;
					temp_LI_pos.WhetherReport = false;
					LI_Positions.push_back(temp_LI_pos);
					EventIndex_Pos[Index_Plus] = LI_Positions_Size;
					EventIndex_Pos[Index_Minus] = LI_Positions_Size;
					LI_Positions_Size++;
					//Index_Minus += 30;
				}
         }
		}
	}
   //cout << "LI: " << LI_Positions.size() << endl;
	static int Count_LI = 0;
	// find LI supporting reads
	
	
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
		//UP_Close_index = Reads[Index].UP_Close.size() - 1;
		if (Reads[Index].Used || !Reads[Index].UP_Far.empty()) continue;
		temp_AbsLoc = Reads[Index].UP_Close[Reads[Index].UP_Close.size() - 1].AbsLoc;
		if (EventIndex_Pos[temp_AbsLoc] == -1) continue;
		Reads[Index].Used = true;
		if (Reads[Index].MatchedD == Plus) {
			LI_Positions[EventIndex_Pos[temp_AbsLoc]].Plus_Reads.push_back(Index);
		}
		else {
			LI_Positions[EventIndex_Pos[temp_AbsLoc]].Minus_Reads.push_back(Index);	
		}
	}
	
	vector <SPLIT_READ> temp_Plus_Reads, temp_Minus_Reads;
	
	bool temp_BalancedPlus_Plus, temp_BalancedPlus_Minus, temp_BalancedMinus_Plus, temp_BalancedMinus_Minus;
	short temp_LengthStr;
	for (unsigned LI_index = 0; LI_index < LI_Positions.size(); LI_index++) {
		if (LI_Positions[LI_index].Minus_Reads.empty() || LI_Positions[LI_index].Plus_Reads.empty()) continue;
		temp_BalancedPlus_Plus = false;
		temp_BalancedPlus_Minus = false;
		temp_BalancedMinus_Plus = false;
		temp_BalancedMinus_Minus = false;
		temp_Plus_Reads.clear();
		temp_Minus_Reads.clear();
		//cout << "Here: " << LI_index << "\t" << LI_Positions[LI_index].Minus_Reads.size() << "\t" << LI_Positions[LI_index].Plus_Reads.size() << endl;
		for (int i = 0; i < LI_Positions[LI_index].Minus_Reads.size(); i++) {
			temp_Minus_Reads.push_back(Reads[LI_Positions[LI_index].Minus_Reads[i]]);
		}
		//CheckConsistancy(temp_Minus_Reads);
		for (int i = 0; i < LI_Positions[LI_index].Minus_Reads.size(); i++) {
			UP_Close_index = temp_Minus_Reads[i].UP_Close.size() - 1;
			temp_LengthStr = temp_Minus_Reads[i].UP_Close[UP_Close_index].LengthStr;
			if ((float)temp_LengthStr > temp_Minus_Reads[i].ReadLength * 0.5)
				temp_BalancedMinus_Plus = true;
			else if ((float)temp_LengthStr < temp_Minus_Reads[i].ReadLength * 0.5)
				temp_BalancedMinus_Minus = true;
		}	
		for (int i = 0; i < LI_Positions[LI_index].Plus_Reads.size(); i++) {
			temp_Plus_Reads.push_back(Reads[LI_Positions[LI_index].Plus_Reads[i]]);
		}
		//CheckConsistancy(temp_Plus_Reads);
		for (int i = 0; i < LI_Positions[LI_index].Plus_Reads.size(); i++) {
			UP_Close_index = temp_Plus_Reads[i].UP_Close.size() - 1;
			temp_LengthStr = temp_Plus_Reads[i].UP_Close[UP_Close_index].LengthStr;
			if ((float)temp_LengthStr > temp_Plus_Reads[i].ReadLength * 0.5)
				temp_BalancedPlus_Plus = true;
			else if ((float)temp_LengthStr < temp_Plus_Reads[i].ReadLength * 0.5)
				temp_BalancedPlus_Minus = true;
		}
		
		unsigned int NumSupportPerTagPlus[VectorTag.size()];
		unsigned int NumSupportPerTagMinus[VectorTag.size()];
		for (short i = 0; i < VectorTag.size(); i++) {
		   NumSupportPerTagPlus[i] = 0;
			NumSupportPerTagMinus[i] = 0;
		}
		for (int i = 0; i < temp_Minus_Reads.size(); i++) {
			for (short j = 0; j < VectorTag.size(); j++) {
			   if (temp_Minus_Reads[i].Tag == VectorTag[j]) {
					NumSupportPerTagMinus[j]++;
					break;
				}	
			}
		}
		for (int i = 0; i < temp_Plus_Reads.size(); i++) {
			for (short j = 0; j < VectorTag.size(); j++) {
			   if (temp_Plus_Reads[i].Tag == VectorTag[j]) {
					NumSupportPerTagPlus[j]++;
					break;
				}	
			}
		}
		bool SupportedByOneSample = false;
		for (short j = 0; j < VectorTag.size(); j++) {
			//cout << NumSupportPerTagPlus[j] << "\t" << NumSupportPerTagMinus[j] << endl;
			if (NumSupportPerTagPlus[j] > 0 && NumSupportPerTagMinus[j] > 0) {
				SupportedByOneSample = true;
				break;
			}
		}
		
		short PositiveBool = 0;
		if (temp_BalancedPlus_Plus) PositiveBool++;
		if (temp_BalancedPlus_Minus) PositiveBool++;
		if (temp_BalancedMinus_Plus) PositiveBool++;
		if (temp_BalancedMinus_Minus) PositiveBool++;
		
		if (SupportedByOneSample && PositiveBool >= 3) 
		{
			//if (LI_Positions[LI_index].WhetherReport) 
			{
				LargeInsertionOutf << "########################################################" << endl;
				LargeInsertionOutf << Count_LI++ << "\tLI\tChrID " << temp_Plus_Reads[0].FragName << "\t" 
				<< LI_Positions[LI_index].Plus_Pos - SpacerBeforeAfter + 1 << "\t" 
				<< temp_Plus_Reads.size() << "\t"
				<< LI_Positions[LI_index].Minus_Pos - SpacerBeforeAfter  + 1 << "\t" 
				<< temp_Minus_Reads.size() << endl;
				
				LargeInsertionOutf <<  ( CurrentChr.substr(LI_Positions[LI_index].Plus_Pos - ReportLength + 1, ReportLength) )
				<< Cap2Low( CurrentChr.substr(LI_Positions[LI_index].Plus_Pos + 1, ReportLength ) ) << endl;
				for (int i = 0; i < temp_Plus_Reads.size(); i++) {
					UP_Close_index = temp_Plus_Reads[i].UP_Close.size() - 1;
					temp_LengthStr = temp_Plus_Reads[i].UP_Close[UP_Close_index].LengthStr;
					for (int j = 0; j < ReportLength - temp_LengthStr; j++) {
						LargeInsertionOutf << " ";
					}
					LargeInsertionOutf << ReverseComplement(temp_Plus_Reads[i].UnmatchedSeq) << endl;
				}
				
				LargeInsertionOutf << "--------------------------------------------------------" << endl;
				//LargeInsertionOutf << "-\t" << minus_LI_Pos[Index_Minus].NumReads << endl;
				LargeInsertionOutf << Cap2Low( CurrentChr.substr(LI_Positions[LI_index].Minus_Pos - ReportLength, ReportLength) )
				<< ( CurrentChr.substr(LI_Positions[LI_index].Minus_Pos, ReportLength ) ) << endl;
				for (int i = 0; i < temp_Minus_Reads.size(); i++) {
					UP_Close_index = temp_Minus_Reads[i].UP_Close.size() - 1;
					temp_LengthStr = temp_Minus_Reads[i].UP_Close[UP_Close_index].LengthStr;
					for (int j = 0; j < ReportLength + temp_LengthStr - temp_Minus_Reads[i].ReadLength; j++) {
						LargeInsertionOutf << " ";
					}
					LargeInsertionOutf << (temp_Minus_Reads[i].UnmatchedSeq) << endl;
				}	
			}
			
		}
		
		//LI_Positions[LI_index].WhetherReport = true; 
	}
	// output
	delete[] plus_LI_Pos ;
	delete[] minus_LI_Pos ;
	delete[] EventIndex_Pos;
	
	cout << "Breakpoints for large insertions (LI): " << Count_LI << "\n\n";
}


void SortOutputRest(const string & CurrentChr, vector <SPLIT_READ> & Reads, vector <SPLIT_READ> & BP_Reads, ofstream & Outf_Rest) {
	SPLIT_READ one_BP_read;
	string HalfMapped, HalfUnmapped;
	int HalfMappedIndex, HalfUnmappedIndex;
	unsigned UP_Close_index;
	unsigned temp_AbsLoc;
	//cout << "1" << endl;
	// find LI combinations
	uint8_t * plus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	uint8_t * minus_LI_Pos = new uint8_t[CurrentChr.size() + 1];
	for (unsigned i = 0; i < CurrentChr.size() + 1; i++) {
		plus_LI_Pos[i] = 0;
		minus_LI_Pos[i] = 0;
	}
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
		if (Reads[Index].Found || Reads[Index].Used  || !Reads[Index].UP_Far.empty()) continue;
		UP_Close_index = Reads[Index].UP_Close.size() - 1;
		temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
		if (plus_LI_Pos[temp_AbsLoc] < NumRead2ReportCutOff_BP) {
			if (Reads[Index].MatchedD == Plus)
				plus_LI_Pos[temp_AbsLoc]++;			
		}
      if (minus_LI_Pos[temp_AbsLoc] < NumRead2ReportCutOff_BP) {
			if (Reads[Index].MatchedD == Minus)
				minus_LI_Pos[temp_AbsLoc]++;	
		}
	}
	vector <Rest_Pos> Rest_Positions;
	Rest_Pos temp_Rest_pos;
	
	//cout << "2" << endl;
	
	bool SkipThisPos;
	for (int Index = SpacerBeforeAfter; Index < CurrentChr.size() - SpacerBeforeAfter; Index++) {
		SkipThisPos = false;
		//for (int MaskIndex = Index + 10; MaskIndex >= Index - 10; MaskIndex--) {
		//	if (CurrentChrMask[MaskIndex] == 'B') {
		//		Index = MaskIndex + 10;
		//		SkipThisPos = true;
		//		break;
		//	}
		//}
		if (SkipThisPos == true) continue;
		if (plus_LI_Pos[Index] >= NumRead2ReportCutOff_BP) {
			temp_Rest_pos.Strand = Plus;
			temp_Rest_pos.Pos = Index;
			Rest_Positions.push_back(temp_Rest_pos);
		}
		if (minus_LI_Pos[Index] >= NumRead2ReportCutOff_BP) {
			temp_Rest_pos.Strand = Minus;
			temp_Rest_pos.Pos = Index;
			Rest_Positions.push_back(temp_Rest_pos);
		}
	}
	
	// find supporting reads
	for (unsigned Index = 0; Index < Reads.size(); Index++) {
		if (Reads[Index].Used || !Reads[Index].UP_Far.empty()) continue;
		UP_Close_index = Reads[Index].UP_Close.size() - 1;
		temp_AbsLoc = Reads[Index].UP_Close[UP_Close_index].AbsLoc;
		for (unsigned Pos_index = 0; Pos_index < Rest_Positions.size(); Pos_index++) {
			if (Reads[Index].MatchedD == Rest_Positions[Pos_index].Strand) {
				if (temp_AbsLoc == Rest_Positions[Pos_index].Pos) {
					Reads[Index].Used = true;
					Rest_Positions[Pos_index].Pos_Reads.push_back(Index);	 // copy index to save memory
				}
			}
		}
	}
	
	//cout << "Other unassigned breakpoints (BP): " << Rest_Positions.size() << "\n\n";
	int Count_BP = 0;
	bool temp_BalancedPlus, temp_BalancedMinus;
	short temp_LengthStr;
	vector <SPLIT_READ> temp_Pos_Reads;
	for (unsigned LI_index = 0; LI_index < Rest_Positions.size(); LI_index++) {
		temp_Pos_Reads.clear();
		temp_BalancedPlus = false;
		temp_BalancedMinus = false;
		for (int i = 0; i < Rest_Positions[LI_index].Pos_Reads.size(); i++) {
			temp_Pos_Reads.push_back(Reads[Rest_Positions[LI_index].Pos_Reads[i]]);
			UP_Close_index = Reads[Rest_Positions[LI_index].Pos_Reads[i]].UP_Close.size() - 1;
			temp_LengthStr = Reads[Rest_Positions[LI_index].Pos_Reads[i]].UP_Close[UP_Close_index].LengthStr;
			if ((float)temp_LengthStr > Reads[Rest_Positions[LI_index].Pos_Reads[i]].ReadLength * 0.5)
				temp_BalancedPlus = true;
			else if ((float)temp_LengthStr < Reads[Rest_Positions[LI_index].Pos_Reads[i]].ReadLength * 0.5)
				temp_BalancedMinus = true;
		}	
		if (temp_BalancedPlus && temp_BalancedMinus) {
			Count_BP++;
			if (Rest_Positions[LI_index].Strand == Plus) {
				Outf_Rest << "########################################################" << endl;
				Outf_Rest << "ChrID " << temp_Pos_Reads[0].FragName << "\t" 
				<< Rest_Positions[LI_index].Pos - SpacerBeforeAfter + 1 << "\t" 
				<< Rest_Positions[LI_index].Pos_Reads.size() << "\t+" << endl;
				
				Outf_Rest <<  ( CurrentChr.substr(Rest_Positions[LI_index].Pos - ReportLength + 1, ReportLength) )
				<< Cap2Low( CurrentChr.substr(Rest_Positions[LI_index].Pos + 1, ReportLength ) ) << endl;
				HalfMappedIndex = 0;
				HalfUnmappedIndex = 0;
				for (int i = 0; i < temp_Pos_Reads.size(); i++) {
					UP_Close_index = temp_Pos_Reads[i].UP_Close.size() - 1;
					temp_LengthStr = temp_Pos_Reads[i].UP_Close[UP_Close_index].LengthStr;
					if (temp_LengthStr > temp_Pos_Reads[HalfMappedIndex].UP_Close[temp_Pos_Reads[HalfMappedIndex].UP_Close.size() - 1].LengthStr)
						HalfMappedIndex = i;
					if (temp_Pos_Reads[i].ReadLength - temp_LengthStr > temp_Pos_Reads[HalfUnmappedIndex].ReadLength - temp_Pos_Reads[HalfUnmappedIndex].UP_Close[temp_Pos_Reads[HalfUnmappedIndex].UP_Close.size() - 1].LengthStr)
						HalfUnmappedIndex = i;
					for (int j = 0; j < ReportLength - temp_LengthStr; j++) {
						Outf_Rest << " ";
					}
					Outf_Rest << ReverseComplement(temp_Pos_Reads[i].UnmatchedSeq)
					<< "\t" << temp_Pos_Reads[i].MatchedD  
					<< "\t" << temp_Pos_Reads[i].MatchedRelPos 
					<< "\t" << temp_Pos_Reads[i].MS
					<< "\t" << temp_Pos_Reads[i].Tag 
					<< "\t" <<  temp_Pos_Reads[i].Name << endl;
					//temp_Pos_Reads[0].
					//BP_Reads.push_back();
				}
				
			}
			else {
				Outf_Rest << "########################################################" << endl;
				Outf_Rest << "ChrID " << temp_Pos_Reads[0].FragName << "\t" 
				<< Rest_Positions[LI_index].Pos - SpacerBeforeAfter + 1 << "\t" 
				<< Rest_Positions[LI_index].Pos_Reads.size() << "\t-" << endl;
				Outf_Rest << Cap2Low( CurrentChr.substr(Rest_Positions[LI_index].Pos - ReportLength, ReportLength) )
				<< ( CurrentChr.substr(Rest_Positions[LI_index].Pos, ReportLength ) ) << endl;
				for (int i = 0; i < temp_Pos_Reads.size(); i++) {
					UP_Close_index = temp_Pos_Reads[i].UP_Close.size() - 1;
					temp_LengthStr = temp_Pos_Reads[i].UP_Close[UP_Close_index].LengthStr;
					for (int j = 0; j < ReportLength + temp_LengthStr - temp_Pos_Reads[i].ReadLength; j++) {
						Outf_Rest << " ";
					}
					Outf_Rest << (temp_Pos_Reads[i].UnmatchedSeq)
					<< "\t" << temp_Pos_Reads[i].MatchedD  
					<< "\t" << temp_Pos_Reads[i].MatchedRelPos 
					<< "\t" << temp_Pos_Reads[i].MS
					<< "\t" << temp_Pos_Reads[i].Tag 
					<< "\t" <<  temp_Pos_Reads[i].Name << endl;				}					 
			}
		} 
	}
	delete[] plus_LI_Pos ;
	delete[] minus_LI_Pos ;
	cout << "Other unassigned breakpoints (BP): " << Count_BP << "\n\n";
}

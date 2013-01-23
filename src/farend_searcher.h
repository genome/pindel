/*
 * This File is part of Pindel; a program to locate genomic variation.
 * https://trac.nbic.nl/pindel/
 *
 *   Copyright (C) 2011 Kai Ye
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FARENDSEARCHER_H_
#define FARENDSEARCHER_H_

class FarEndSearcher {
public:
	FarEndSearcher(const std::string & CurrentChr_in, SPLIT_READ & Temp_One_Read_in);
    void GetFarEnd(const int &in_start, const int &in_end);
    void GetFarEnd_OtherStrand(const short &RangeIndex);
    void GetFarEnd_SingleStrandUpStream(const short &RangeIndex);
    void GetFarEnd_SingleStrandDownStream(const short &RangeIndex);
    //void SearchFarEndAtPos( const std::string& chromosome, SPLIT_READ& read, const unsigned int searchCenter, const unsigned int range);
	virtual ~FarEndSearcher();

private:
	void GetFarEnd_General(const int &in_start, const int &in_end, const bool &UseRangeIndex, const short &RangeIndex);
	const std::string* CurrentChr;
	SPLIT_READ* Temp_One_Read;
};

void SearchFarEndAtPos( const std::string& chromosome, SPLIT_READ& Temp_One_Read, const unsigned int SearchCenter, const unsigned int Range);
void SearchFarEndAtPos( const std::string& chromosome, SPLIT_READ& Temp_One_Read, const std::vector <SearchWindow> & Regions );

#endif /* FARENDSEARCHER_H_ */

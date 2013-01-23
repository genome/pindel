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

#ifndef REFREADER_H
#define	REFREADER_H

class RefReader
{
public:
  void ReadChr(const std::string & ChrName) ;
  void ReadSeq(bool WhetherBuildUp) ;
  RefReader(std::ifstream & inf_Seq_in, std::string & TheInput_in) ;
  ~RefReader() ;

private:
  RefReader(const RefReader&);

  std::ifstream* inf_Seq;
  std::string* TheInput;
  void CopyThisSeq() ;
  void SkipThisSeq() ;
};


#endif	/* REFREADER_H */

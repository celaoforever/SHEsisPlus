/*
 * CreatHtmlTable.cpp
 *
 *  Created on: Sep 21, 2014
 *      Author: ionadmin
 */

#include "CreatHtmlTable.h"
namespace SHEsis{
CreatHtmlTable::CreatHtmlTable() {
	// TODO Auto-generated constructor stub

}

CreatHtmlTable::~CreatHtmlTable() {
	// TODO Auto-generated destructor stub
}

void CreatHtmlTable::createTable(std::string id){
	res<<"<table id=\""<<id<<"\">\n";
};
void CreatHtmlTable::createTable(){
	res<<"<table>\n";
};
std::string CreatHtmlTable::getTable(){
	res<<"</table>\n";
	return res.str();
}
void CreatHtmlTable::addHeadRow(std::vector<std::string>& v ){
	res<<"<tr>\n";
	for(int i=0;i<v.size();i++){
		res<<"<th>"<<v[i]<<"</th>\n";
	}
	res<<"</tr>\n";
}
void CreatHtmlTable::addDataRow(std::vector<std::string>& v){
	res<<"<tr>\n";
	for(int i=0;i<v.size();i++){
		res<<"<td>"<<v[i]<<"</td>\n";
	}
	res<<"</tr>\n";
}
};


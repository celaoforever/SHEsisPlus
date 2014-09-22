/*
 * CreatHtmlTable.h
 *
 *  Created on: Sep 21, 2014
 *      Author: ionadmin
 */

#ifndef CREATHTMLTABLE_H_
#define CREATHTMLTABLE_H_
#include <string>
#include <sstream>
#include <vector>
namespace SHEsis {
class CreatHtmlTable {
 public:
  CreatHtmlTable();
  virtual ~CreatHtmlTable();
  void createTable(std::string id);
  void createTable();
  std::string getTable();
  void addHeadRow(std::vector<std::string>& v);
  void addDataRow(std::vector<std::string>& v);

 private:
  std::string id;
  std::stringstream res;
};
}
#endif /* CREATHTMLTABLE_H_ */

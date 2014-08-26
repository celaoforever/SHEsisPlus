/*
 * IndexingVariables.cpp
 *
 *  Created on: Aug 18, 2014
 *      Author: ada
 */


#include "IndexingVariables.h"


namespace SHEsis{
	void IndexingVariables::reset(){
		if(id.size()>0)
			this->id.clear();
		if(enumeration.size()>0)
			this->enumeration.clear();
	}


	void IndexingVariables::add(std::string key, boost::shared_ptr<int[]> sizes)
	{
	    int size = GeneralIndexingGetSize(sizes);
	    //std::cout<<key<<"\n";
	    //BOOST_ASSERT(size);
	    boost::shared_ptr<int[]> variables(new int[size+1]);
	    variables[0]=size;
	    for(int i = 1; i <= size; i++)
	        variables[i] = offset + i;
	    offset += size;
	    boost::shared_ptr<ArrayStorage> spAS(new ArrayStorage(sizes, variables));
	    hm[key]=spAS;
	};
	void IndexingVariables::printHmKey(){
		boost::unordered_map<std::string,boost::shared_ptr<ArrayStorage> >::iterator iter;
		std::cout<<"\nkey of hm:\n";
		for(iter=this->hm.begin();iter!=this->hm.end();iter++)
		{
			std::cout<<iter->first<<":";
			for(int i=1;i<=iter->second->sizes[0];i++){
				std::cout<<iter->second->sizes[i]<<",";
			}
			std::cout<<"\n";
		}
	}
	int IndexingVariables::getEnumeration(int ident){
		BOOST_ASSERT(ident>=1);
		return this->id[ident];
	};

	int IndexingVariables::getEnumeration(std::string key, boost::shared_ptr<int[]> indices){
		BOOST_ASSERT(hm.find(key) != hm.end());
		int ident=SHEsisABS((this->hm[key])->get(indices));
		int size=0;
		if(id.end() == id.find(ident)){
			size = id.size()+1;
			id[ident]=size;
			enumeration[size]=ident;
			return size;
		}else{
			return id[ident];
		}
	};

	std::string IndexingVariables::getEnumerationILP( std::string key, boost::shared_ptr<int[]> indices){
		int ident=SHEsisABS((this->hm[key])->get(indices));
		int size=0;
		std::stringstream res;
		if(id.end() == id.find(ident)){
			size=id.size()+1;
			id[ident]=size;
			enumeration[size]=ident;
			res<<"ilp."<<size;
			return res.str();
		}else
		{
			res<<"ilp."<<id[ident];
			return res.str();
		}
	};

	int IndexingVariables::getId(int enumerat){
		BOOST_ASSERT(enumerat>=1);
		return enumeration[enumerat];
	}

	int IndexingVariables::getEvalutatedId(std::string key, boost::shared_ptr<int[]> indices){
		return this->hm[key]->get(indices);
	};

	double IndexingVariables::getEvalutatedIdILP(std::string key, boost::shared_ptr<int[]> indices){
		int ident=this->hm[key]->get(indices);
		if(id.end() != id.find(ident)){
			std::stringstream k;
			k<<"ilp."<<id[ident];
			return ILPsolution[k.str()];
		}else
		{
			return NULL;
		}
	};

	void IndexingVariables::set( std::string key, boost::shared_ptr<int[]> variables, boost::shared_ptr<int[]> sizes ){
		boost::shared_ptr<ArrayStorage> spAS(new ArrayStorage(sizes, variables));
		hm[key]=spAS;
	};

	void IndexingVariables::set(boost::unordered_map<std::string,boost::shared_ptr<ArrayStorage> > as, int offset){
		hm.clear();
		hm.insert(as.begin(),as.end());
		this->offset=offset;
	}

	boost::shared_ptr<int[]> IndexingVariables::getObjectSizes(std::string key){
		if(hm.end() != hm.find(key))
			return hm[key]->getSizes();
		boost::shared_ptr<int[]> res(new int[1]);
		res[0]=-1;
		return res;
	};

	int IndexingVariables::size()
	{
		return hm.size();
	};

	int IndexingVariables::numberOfVariables()
	{
		return this->offset-1;
	}

	void IndexingVariables::setParities(boost::multi_array<int,1> parities, boost::multi_array<int, 1> outsider){
		//TO-DO
	}

	void IndexingVariables::setValueILP(std::string key, double value){
		ILPsolution[key]=value;
	};


	void IndexingVariables::show(){};
	int IndexingVariables::test(){return id.size();};

	int IndexingVariables::getEnumerationSize(){
		return enumeration.size();
	};

	bool IndexingVariables::containsKey(std::string key){
		return (hm.end() != hm.find(key));
	};

	int IndexingVariables::getIdSize(){
		return id.size();
	};

}


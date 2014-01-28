/*
 * DataPointDefinition.h
 *
 *  Created on: Sep 24, 2012
 *      Author: aspataru
 */

#ifndef DATAPOINTDEFINITION_H_
#define DATAPOINTDEFINITION_H_

#include "EventFilter/Utilities/interface/JsonMonitorable.h"
#include "EventFilter/Utilities/interface/JsonSerializable.h"
#include <string>
#include <vector>

namespace jsoncollector {

class JsonMonConfig;

class DataPointDefinition: public JsonSerializable {

public:
	DataPointDefinition() {}

	//DataPointDefinition(const std::vector<std::string>& names, const std::vector<std::string>& operations);

	virtual ~DataPointDefinition() {}

	/**
	 * JSON serialization procedure for this class
	 */
	virtual void serialize(Json::Value& root) const;
	/**
	 * JSON deserialization procedure for this class
	 */
	virtual void deserialize(Json::Value& root);
	/**
	 * Returns true if the legend_ has elements
	 */
	bool isPopulated() const;
	/**
	 * Returns a LegendItem object ref at the specified index
	 */
	std::vector<std::string> const& getNames() {return varNames_;}
	std::vector<std::string> const& getOperations() {return opNames_;}

	/**
	 * Loads a DataPointDefinition from a specified reference
	 */
	static bool getDataPointDefinitionFor(std::string& defFilePath, DataPointDefinition& dpd);

	OperationType getOperationFor(unsigned int index);

	std::string & getDefFilePath() {return defFilePath_;}
	//void populateMonConfig(std::vector<JsonMonConfig>& monConfig);

	//known JSON operation names
	static const std::string SUM;
	static const std::string AVG;
	static const std::string SAME;
	static const std::string HISTO;
	static const std::string CAT;

	// JSON field names
	static const std::string LEGEND;
	static const std::string PARAM_NAME;
	static const std::string OPERATION;

private:
	std::vector<std::string> varNames_;
	std::vector<std::string> opNames_;
	std::string defFilePath_;
};
}

#endif /* DATAPOINTDEFINITION_H_ */

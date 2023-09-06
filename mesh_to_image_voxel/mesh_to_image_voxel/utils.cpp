
#include <intrin.h>  
#include<fstream>
#include<stdio.h>
#include<math.h>
#include<filesystem>
#include<iostream>
namespace fs = std::experimental::filesystem;

void ergodicfile(const std::string &root, std::vector<std::string> &vecFile, const std::string &fileType)
{
	for (auto &dir : fs::directory_iterator(root))
	{
		if (fs::is_directory(dir))
		{
			auto path = dir.path().string();
			ergodicfile(path, vecFile, fileType);
		}
		else
		{
			auto file = dir.path().string();
			if (file.find(fileType) != -1)
			{
				vecFile.push_back(file);
			}
		}
	}
}
//This is the customized trie data structure to return the previous word id if the search exist.
#pragma once
#ifndef CTRIE_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define CTRIE_H
#include "stdafx.h"
class CTrieNode
{
public:
    std::unordered_map<int, CTrieNode*> children;
    bool isLeaf;
    int firstID;
    CTrieNode()
    {
        firstID = -1;
        isLeaf = false;
    }
};

class CTrie
{
public:
    CTrieNode* root;
    CTrie()
    {
        root = new CTrieNode();
    }

    void insert(std::vector<int> nbrSet, int ID);
    bool search(std::vector<int> nbrSet, int& firstID);
};

void CTrie::insert(std::vector<int> nbrSet, int ID)
{
    CTrieNode* current = root;

    for (int i = 0; i < nbrSet.size(); ++i)
    {
        int nbr = nbrSet[i];
        if (current->children.find(nbr) == current->children.end())
        {
            current->children.insert(make_pair(nbr, new CTrieNode()));
        }
        current = current->children[nbr];
    }
    current->firstID = ID;
    current->isLeaf = true;
}

bool CTrie::search(std::vector<int> nbrSet, int& firstID)
{
    CTrieNode* current = root;
    for (int i = 0; i < nbrSet.size(); ++i)
    {
        int nbr = nbrSet[i];
        if (current->children.find(nbr) != current->children.end())
        {
            current = current->children[nbr];
        }
        else
        {
            return false;
        }
    }
    firstID = current->firstID;
    return current->isLeaf;
}
#endif


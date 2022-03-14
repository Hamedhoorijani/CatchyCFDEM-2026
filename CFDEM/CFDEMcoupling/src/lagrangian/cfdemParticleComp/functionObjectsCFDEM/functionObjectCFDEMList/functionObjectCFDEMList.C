/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "functionObjectCFDEMList.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "argList.H"
#include "timeControlFunctionObject.H"
#include "dictionaryEntry.H"
#include "stringOps.H"
#include "Tuple2.H"
#include "etcFiles.H"
#include "IOdictionary.H"
#include "wordAndDictionary.H"


/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

Foam::fileName Foam::functionObjectCFDEMList::functionObjectDictPath
(
    "caseDicts/postProcessing"
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::functionObjectCFDEM* Foam::functionObjectCFDEMList::remove
(
    const word& key,
    label& oldIndex
)
{
    functionObjectCFDEM* ptr = 0;

    // Find index of existing functionObject
    HashTable<label>::iterator fnd = indices_.find(key);

    if (fnd != indices_.end())
    {
        oldIndex = fnd();

        // Retrieve the pointer and remove it from the old list
        ptr = this->set(oldIndex, 0).ptr();
        indices_.erase(fnd);
    }
    else
    {
        oldIndex = -1;
    }

    return ptr;
}


void Foam::functionObjectCFDEMList::listDir
(
    const fileName& dir,
    HashSet<word>& foMap
)
{
    // Search specified directory for functionObject configuration files
    {
        fileNameList foFiles(fileHandler().readDir(dir));
        forAll(foFiles, f)
        {
            if (foFiles[f].ext().empty())
            {
                foMap.insert(foFiles[f]);
            }
        }
    }

    // Recurse into sub-directories
    {
        fileNameList foDirs(fileHandler().readDir(dir, fileType::directory));
        forAll(foDirs, fd)
        {
            listDir(dir/foDirs[fd], foMap);
        }
    }
}


void Foam::functionObjectCFDEMList::list()
{
    HashSet<word> foMap;

    fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

    forAll(etcDirs, ed)
    {
        listDir(etcDirs[ed], foMap);
    }

    Info<< nl
        << "Available configured functionObjects:"
        << foMap.sortedToc()
        << nl;
}


Foam::fileName Foam::functionObjectCFDEMList::findDict
(
    const word& funcName,
    const word& region
)
{
    // First check if there is a functionObject dictionary file in the
    // region system directory
    {
        const fileName dictFile
        (
            stringOps::expand("$FOAM_CASE")/"system"/region/funcName
        );

        if (isFile(dictFile))
        {
            return dictFile;
        }
    }

    // Next, if the region is specified, check if there is a functionObject
    // dictionary file in the global system directory
    if (region != word::null)
    {
        const fileName dictFile
        (
            stringOps::expand("$FOAM_CASE")/"system"/funcName
        );

        if (isFile(dictFile))
        {
            return dictFile;
        }
    }

    // Finally, check etc directories
    {
        const fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

        forAll(etcDirs, i)
        {
            const fileName dictFile(search(funcName, etcDirs[i]));

            if (!dictFile.empty())
            {
                return dictFile;
            }
        }
    }

    return fileName::null;
}


void Foam::functionObjectCFDEMList::checkUnsetEntries
(
    const string& funcCall,
    const dictionary& funcArgsDict,
    const dictionary& funcDict,
    const string& context
)
{
    const wordRe unset("<.*>");
    unset.compile();

    forAllConstIter(IDLList<entry>, funcArgsDict, iter)
    {
        if (iter().isStream())
        {
            ITstream& tokens = iter().stream();

            forAll(tokens, i)
            {
                if (tokens[i].isWord())
                {
                    if (unset.match(tokens[i].wordToken()))
                    {
                        FatalIOErrorInFunction(funcDict)
                            << "Essential value for keyword '"
                            << iter().keyword()
                            << "' not set in function entry" << nl
                            << "    " << funcCall.c_str() << nl
                            << "    in " << context.c_str() << nl
                            << "    Placeholder value is "
                            << tokens[i].wordToken()
                            << exit(FatalIOError);
                    }
                }
            }
        }
        else
        {
            checkUnsetEntries(funcCall, iter().dict(), funcDict, context);
        }
    }
}


bool Foam::functionObjectCFDEMList::readFunctionObject
(
    const string& funcCall,
    dictionary& functionsDict,
    const string& context,
    HashSet<word>& requiredFields,
    const word& region
)
{
    word funcName(funcCall);

    int argLevel = 0;
    wordReList args;

    List<Tuple2<word, string>> namedArgs;
    bool namedArg = false;
    word argName;

    word::size_type start = 0;
    word::size_type i = 0;

    for
    (
        word::const_iterator iter = funcCall.begin();
        iter != funcCall.end();
        ++iter
    )
    {
        char c = *iter;

        if (c == '(')
        {
            if (argLevel == 0)
            {
                funcName = funcCall(start, i - start);
                start = i+1;
            }
            ++argLevel;
        }
        else if (c == ',' || c == ')')
        {
            if (argLevel == 1)
            {
                if (namedArg)
                {
                    namedArgs.append
                    (
                        Tuple2<word, string>
                        (
                            argName,
                            funcCall(start, i - start)
                        )
                    );
                    namedArg = false;
                }
                else
                {
                    args.append(wordRe(funcCall(start, i - start)));
                }
                start = i+1;
            }

            if (c == ')')
            {
                if (argLevel == 1)
                {
                    break;
                }
                --argLevel;
            }
        }
        else if (c == '=')
        {
            argName = string::validate<word>(funcCall(start, i - start));
            start = i+1;
            namedArg = true;
        }

        ++i;
    }

    // Strip whitespace from the function name
    string::stripInvalid<word>(funcName);

    // Search for the functionObject dictionary
    fileName path = findDict(funcName, region);

    if (path == fileName::null)
    {
        WarningInFunction
            << "Cannot find functionObject file " << funcName << endl;
        return false;
    }

    // Read the functionObject dictionary
    // IFstream fileStream(path);
    autoPtr<ISstream> fileStreamPtr(fileHandler().NewIFstream(path));
    ISstream& fileStream = fileStreamPtr();

    // Delay processing the functionEntries
    // until after the function argument entries have been added
    entry::disableFunctionEntries = true;
    dictionary funcsDict(funcName, functionsDict, fileStream);
    entry::disableFunctionEntries = false;

    dictionary* funcDictPtr = &funcsDict;

    if (funcsDict.found(funcName) && funcsDict.isDict(funcName))
    {
        funcDictPtr = &funcsDict.subDict(funcName);
    }

    dictionary& funcDict = *funcDictPtr;

    // Store the funcDict as read for error reporting context
    const dictionary funcDict0(funcDict);

    // Insert the 'field' and/or 'fields' and 'objects' entries corresponding
    // to both the arguments and the named arguments
    DynamicList<wordAndDictionary> fieldArgs;
    forAll(args, i)
    {
        fieldArgs.append(wordAndDictionary(args[i], dictionary::null));
    }
    forAll(namedArgs, i)
    {
        if (namedArgs[i].first() == "field")
        {
            IStringStream iss(namedArgs[i].second());
            fieldArgs.append(wordAndDictionary(iss));
        }
        if
        (
            namedArgs[i].first() == "fields"
         || namedArgs[i].first() == "objects"
        )
        {
            IStringStream iss(namedArgs[i].second());
            fieldArgs.append(List<wordAndDictionary>(iss));
        }
    }
    if (fieldArgs.size() == 1)
    {
        funcDict.set("field", fieldArgs[0].first());
        funcDict.merge(fieldArgs[0].second());
    }
    if (fieldArgs.size() >= 1)
    {
        funcDict.set("fields", fieldArgs);
        funcDict.set("objects", fieldArgs);
    }

    // Insert non-field arguments
    forAll(namedArgs, i)
    {
        if
        (
            namedArgs[i].first() != "field"
         && namedArgs[i].first() != "fields"
         && namedArgs[i].first() != "objects"
        )
        {
            IStringStream entryStream
            (
                namedArgs[i].first() + ' ' + namedArgs[i].second() + ';'
            );
            funcDict.set(entry::New(entryStream).ptr());
        }
    }

    // Insert the region name if specified
    if (region != word::null)
    {
        funcDict.set("region", region);
    }

    const word funcCallKeyword = string::validate<word>(funcCall);
    dictionary funcArgsDict;
    funcArgsDict.add(funcCallKeyword, funcDict);

    // Re-parse the funcDict to execute the functionEntries
    // now that the function argument entries have been added
    {
        OStringStream os;
        funcArgsDict.write(os);
        funcArgsDict = dictionary
        (
            funcName,
            functionsDict,
            IStringStream(os.str())()
        );
    }

    // Check for anything in the configuration that has not been set
    checkUnsetEntries(funcCall, funcArgsDict, funcDict0, context);

    // Lookup the field, fields and objects entries from the now expanded
    // funcDict and insert into the requiredFields
    dictionary& expandedFuncDict = funcArgsDict.subDict(funcCallKeyword);
    if (functionObjectCFDEM::debug)
    {
        InfoInFunction
            << nl << incrIndent << indent
            << funcCall << expandedFuncDict
            << decrIndent << endl;
    }
    if (expandedFuncDict.found("field"))
    {
        requiredFields.insert(word(expandedFuncDict.lookup("field")));
    }
    if (expandedFuncDict.found("fields"))
    {
        List<wordAndDictionary> fields(expandedFuncDict.lookup("fields"));
        forAll(fields, i)
        {
            requiredFields.insert(fields[i].first());
        }
    }
    if (expandedFuncDict.found("objects"))
    {
        List<wordAndDictionary> objects(expandedFuncDict.lookup("objects"));
        forAll(objects, i)
        {
            requiredFields.insert(objects[i].first());
        }
    }

    // Merge this functionObject dictionary into functionsDict
    functionsDict.merge(funcArgsDict);
    functionsDict.subDict(funcCallKeyword).name() = funcDict.name();

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectCFDEMList::functionObjectCFDEMList
(
    const Time& t,
    cfdemCloud& sm,
    const bool execution
)
:
    PtrList<functionObjectCFDEM>(),
    digests_(),
    indices_(),
    time_(t),
    parentDict_(t.controlDict()),
    execution_(execution),
    updated_(false),
    particleCloud_(sm)
{}


Foam::functionObjectCFDEMList::functionObjectCFDEMList
(
    const Time& t,
    const dictionary& parentDict,
    cfdemCloud& sm,
    const bool execution
)
:
    PtrList<functionObjectCFDEM>(),
    digests_(),
    indices_(),
    time_(t),
    parentDict_(parentDict),
    execution_(execution),
    updated_(false),
    particleCloud_(sm)
{}


Foam::autoPtr<Foam::functionObjectCFDEMList> Foam::functionObjectCFDEMList::New
(
    const argList& args,
    const Time& runTime,
    dictionary& controlDict,
    cfdemCloud& sm,
    HashSet<word>& requiredFields
)
{
    autoPtr<functionObjectCFDEMList> functionsPtr;

    controlDict.add
    (
        dictionaryEntry("functionsCFDEM", controlDict, dictionary::null)
    );

    dictionary& functionsDict = controlDict.subDict("functionsCFDEM");

    word region = word::null;

    // Set the region name if specified
    if (args.optionFound("region"))
    {
        region = args["region"];
    }

    if
    (
        args.optionFound("dict")
     || args.optionFound("func")
     || args.optionFound("funcs")
    )
    {
        if (args.optionFound("dict"))
        {
            controlDict.merge
            (
                IOdictionary
                (
                    IOobject
                    (
                        args["dict"],
                        runTime,
                        IOobject::MUST_READ_IF_MODIFIED
                    )
                )
            );
        }

        if (args.optionFound("func"))
        {
            readFunctionObject
            (
                args["func"],
                functionsDict,
                "command line " + args.commandLine(),
                requiredFields,
                region
            );
        }

        if (args.optionFound("funcs"))
        {
            wordList funcs(args.optionLookup("funcs")());

            forAll(funcs, i)
            {
                readFunctionObject
                (
                    funcs[i],
                    functionsDict,
                    "command line " + args.commandLine(),
                    requiredFields,
                    region
                );
            }
        }

        functionsPtr.reset(new functionObjectCFDEMList(runTime, controlDict, sm));
    }
    else
    {
        functionsPtr.reset(new functionObjectCFDEMList(runTime, sm));
    }

    functionsPtr->read();

    return functionsPtr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjectCFDEMList::~functionObjectCFDEMList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjectCFDEMList::clear()
{
    PtrList<functionObjectCFDEM>::clear();
    digests_.clear();
    indices_.clear();
    updated_ = false;
}


Foam::label Foam::functionObjectCFDEMList::findObjectID(const word& name) const
{
    forAll(*this, oi)
    {
        if (operator[](oi).name() == name)
        {
            return oi;
        }
    }

    return -1;
}


void Foam::functionObjectCFDEMList::on()
{
    execution_ = true;
}


void Foam::functionObjectCFDEMList::off()
{
    // For safety, also force a read() when execution is turned back on
    updated_ = execution_ = false;
}


bool Foam::functionObjectCFDEMList::status() const
{
    return execution_;
}


bool Foam::functionObjectCFDEMList::start()
{
    bool ok = read();

    if (execution_)
    {
        forAll(*this, oi)
        {
            if (operator[](oi).executeAtStart())
            {
                ok = operator[](oi).execute() && ok;
                ok = operator[](oi).write() && ok;
            }
        }
    }

    return ok;
}


bool Foam::functionObjectCFDEMList::execute()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, oi)
        {
            ok = operator[](oi).execute() && ok;
            ok = operator[](oi).write() && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectCFDEMList::end()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, oi)
        {
            ok = operator[](oi).end() && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectCFDEMList::setTimeStep()
{
    bool set = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        wordList names;

        forAll(*this, oi)
        {
            if (operator[](oi).setTimeStep())
            {
                names.append(operator[](oi).name());
                set = true;
            }
        }

        if (names.size() > 1)
        {
            WarningInFunction << "Multiple function objects (" << names[0];
            for (label i = 1; i < names.size(); ++ i)
            {
                WarningInFunction << ", " << names[i];
            }
            WarningInFunction << ") are setting the time step." << endl;
        }
    }

    return set;
}


Foam::scalar Foam::functionObjectCFDEMList::timeToNextWrite()
{
    scalar result = vGreat;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, oi)
        {
            result = min(result, operator[](oi).timeToNextWrite());
        }
    }

    return result;
}


bool Foam::functionObjectCFDEMList::read()
{
    bool ok = true;
    updated_ = execution_;

    // Avoid reading/initializing if execution is off
    if (!execution_)
    {
        return true;
    }

    // Update existing and add new functionObjects
    const entry* entryPtr = parentDict_.lookupEntryPtr
    (
        "functions",
        false,
        false
    );

    if (entryPtr)
    {
        PtrList<functionObjectCFDEM> newPtrs;
        List<SHA1Digest> newDigs;
        HashTable<label> newIndices;

        label nFunc = 0;

        if (!entryPtr->isDict())
        {
            FatalIOErrorInFunction(parentDict_)
                << "'functions' entry is not a dictionary"
                << exit(FatalIOError);
        }

        const dictionary& functionsDict = entryPtr->dict();

        libs.open
        (
            functionsDict,
            "libs",
            functionObjectCFDEM::dictionaryConstructorTablePtr_
        );

        newPtrs.setSize(functionsDict.size());
        newDigs.setSize(functionsDict.size());

        forAllConstIter(dictionary, functionsDict, iter)
        {
            const word& key = iter().keyword();

            if (!iter().isDict())
            {
                if (key != "libs")
                {
                    IOWarningInFunction(parentDict_)
                        << "Entry " << key << " is not a dictionary" << endl;
                }

                continue;
            }

            const dictionary& dict = iter().dict();
            bool enabled = dict.lookupOrDefault("enabled", true);

            newDigs[nFunc] = dict.digest();

            label oldIndex;
            functionObjectCFDEM* objPtr = remove(key, oldIndex);

            if (objPtr)
            {
                if (enabled)
                {
                    // Dictionary changed for an existing functionObject
                    if (newDigs[nFunc] != digests_[oldIndex])
                    {
                        ok = objPtr->read(dict) && ok;
                    }
                }
                else
                {
                    // Delete the disabled functionObject
                    delete objPtr;
                    objPtr = nullptr;
                    continue;
                }
            }
            else if (enabled)
            {
                autoPtr<functionObjectCFDEM> foPtr;

                FatalError.throwExceptions();
                FatalIOError.throwExceptions();
                try
                {
                    if
                    (
                        dict.found("writeControl")
                     || dict.found("outputControl")
                    )
                    {
                        foPtr.set
                        (
                            new functionObjectsCFDEM::timeControl(key, time_, dict, particleCloud_)
                        );
                    }
                    else
                    {
                        foPtr = functionObjectCFDEM::New(key, time_, dict, particleCloud_);
                    }
                }
                catch (Foam::IOerror& ioErr)
                {
                    Info<< ioErr << nl << endl;
                    ::exit(1);
                }
                catch (Foam::error& err)
                {
                    WarningInFunction
                        << "Caught FatalError " << err << nl << endl;
                }
                FatalError.dontThrowExceptions();
                FatalIOError.dontThrowExceptions();

                if (foPtr.valid())
                {
                    objPtr = foPtr.ptr();
                }
                else
                {
                    ok = false;
                }
            }

            // Insert active functionObjects into the list
            if (objPtr)
            {
                newPtrs.set(nFunc, objPtr);
                newIndices.insert(key, nFunc);
                nFunc++;
            }
        }

        newPtrs.setSize(nFunc);
        newDigs.setSize(nFunc);

        // Updating the PtrList of functionObjects deletes any
        // existing unused functionObjects
        PtrList<functionObjectCFDEM>::transfer(newPtrs);
        digests_.transfer(newDigs);
        indices_.transfer(newIndices);
    }
    else
    {
        PtrList<functionObjectCFDEM>::clear();
        digests_.clear();
        indices_.clear();
    }

    return ok;
}


void Foam::functionObjectCFDEMList::updateMesh(const mapPolyMesh& mpm)
{
    if (execution_)
    {
        forAll(*this, oi)
        {
            operator[](oi).updateMesh(mpm);
        }
    }
}


void Foam::functionObjectCFDEMList::movePoints(const polyMesh& mesh)
{
    if (execution_)
    {
        forAll(*this, oi)
        {
            operator[](oi).movePoints(mesh);
        }
    }
}


// ************************************************************************* //

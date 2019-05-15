#include "llvm/Support/CommandLine.h"

#include "clang/Frontend/CompilerInstance.h"
#include "clang/Frontend/FrontendAction.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Rewrite/Core/Rewriter.h"
#include "clang/Tooling/CommonOptionsParser.h"
#include "clang/Tooling/Tooling.h"

using namespace clang;
using namespace clang::tooling;
using namespace llvm;

static cl::OptionCategory MyToolCategory("My tool options");
static cl::opt<std::string> MyToolOutput("o", cl::desc("Output file"), cl::cat(MyToolCategory));
static cl::opt<bool> LoopsVectorization("lv", cl::desc("Loops vectorization"), cl::cat(MyToolCategory));


const char loops[] = "for (Int _j = 0; _j < _size.y; _j++) {\n"
                     "for (Int _i = 0; _i < _size.x; _i++) {\n"
                     "_evaluate(_i, _j, Real());\n"
                     "}\n"
                     "}\n";


const char vecLoops[] =
        "for (Int _j = 0; _j < _size.y; _j++) {\n"
        "Int _i = 0;\n"
        "for (; _i < _size.x / 4 * 4; _i+=4) {\n"
        "_evaluate(_i, _j, Real4());\n"
        "}\n"
        "for (; _i < _size.x; _i++) {\n"
        "_evaluate(_i, _j, Real());\n"
        "}\n"
        "}\n";

std::string escapedName(std::string name) {
    for (char& c : name) {
        if (!std::isalnum(c)) {
            switch (c) {
            case '-': c = 'm'; break;
            case '+': c = 'p'; break;
            default: c = '_';
            }
        }
    }
    return '_' + name;
}

std::string temporaryName(const std::string& name) {
    return 't' + escapedName(name);
}

std::string getText(Stmt* stmt, Rewriter& rewriter) {
    assert(stmt);
    return Lexer::getSourceText(
            CharSourceRange::getTokenRange(stmt->getSourceRange()),
            rewriter.getSourceMgr(),
            LangOptions());
}

class SequenceVisitor : public RecursiveASTVisitor<SequenceVisitor> {
public:
    SequenceVisitor(Rewriter& rewriter)
        : mRewriter(rewriter)
    {

    }

    bool TraverseStmt(Stmt* e) {
        if (e && isa<Expr>(e) &&
                (cast<Expr>(e)->getType().getAsString() == "struct Var" ||
                 cast<Expr>(e)->getType().getAsString() == "struct Array"))
        {
            std::string text = getText(e, mRewriter);
            
            mRewriter.ReplaceText(e->getSourceRange(), escapedName(text));

            rhsVars.insert(text);
            
            return true;
        } else if (e && isa<Expr>(e) && cast<Expr>(e)->getType().getAsString() == "struct Temp") {
            std::string tempName = getText(e, mRewriter);
            temps.insert(tempName);
            
            mRewriter.ReplaceText(e->getSourceRange(), escapedName(tempName));

            return true;
        } else if (e && isa<CXXOperatorCallExpr>(e) &&
                   cast<CXXOperatorCallExpr>(e)->getType().getAsString() == "struct Assign") {
            auto ee = cast<CXXOperatorCallExpr>(e);
            assert(ee->getNumArgs() == 2);

            Expr* lhs = ee->getArg(0);
        
            assert(lhs->getType().getAsString() == "struct Array" ||
                lhs->getType().getAsString() == "struct Var" ||
                lhs->getType().getAsString() == "struct Temp");

            std::string lhsText = getText(lhs, mRewriter);
        
            mRewriter.ReplaceText(lhs->getSourceRange(), escapedName(lhsText));

            if (lhs->getType().getAsString() == "struct Array" || 
                lhs->getType().getAsString() == "struct Var") 
            {
                lhsVars.insert(lhsText);
            }
            
            return RecursiveASTVisitor<SequenceVisitor>::TraverseStmt(ee->getArg(1));
        } else {
            return RecursiveASTVisitor<SequenceVisitor>::TraverseStmt(e);
        }

    }

    bool VisitDeclRefExpr(DeclRefExpr* e) {
        if (e->getNameInfo().getAsString() == "operator,") {
            mRewriter.ReplaceText(e->getSourceRange(), ";");
        }

        return true;
    }
    
    std::set<std::string> allVars() {
        std::set<std::string> res;
        res.insert(lhsVars.begin(), lhsVars.end());
        res.insert(rhsVars.begin(), rhsVars.end());
        return res;
    }

    std::set<std::string> rhsVars; // each use of Var
    std::set<std::string> lhsVars; // each use of Var
    std::set<std::string> temps; // each use of Temp

private:
    Rewriter& mRewriter;
};


class MyASTVisitor : public RecursiveASTVisitor<MyASTVisitor> {
public:
    MyASTVisitor(Rewriter& rewriter)
        : mRewriter(rewriter)
    {
    }


    bool VisitVarDecl(VarDecl *v) {
        if (v->getType().getAsString() == "struct Temp") {
//            auto r = v->getTypeSourceInfo()->getTypeLoc().getSourceRange();
//            mRewriter.ReplaceText(r, "Real");
            mRewriter.RemoveText(v->getSourceRange());
        }
        return true;
    }

    bool VisitFunctionDecl(FunctionDecl *f) {
        bool exprDetected = false;
        if (f->getReturnType().getAsString() == "struct Expr") {
            mRewriter.ReplaceText(f->getReturnTypeSourceRange(), "auto");
        }

        int templatesNum = 0;

        for (unsigned i = 0; i < f->getNumParams(); i++) {
            ParmVarDecl* p = f->getParamDecl(i);
            if (p->getType().getAsString() == "struct Expr") {
                mRewriter.ReplaceText(p->getTypeSourceInfo()->getTypeLoc().getSourceRange(),
                                      "E" + std::to_string(templatesNum));
                templatesNum += 1;
            }
        }
        if (templatesNum > 0) {
            mRewriter.InsertText(f->getBeginLoc(), "template <", true, true);
            for (int i = 0; i < templatesNum; i++) {
                mRewriter.InsertText(f->getBeginLoc(), "typename E" + std::to_string(i), true, true);
                if (i != templatesNum - 1) {
                    mRewriter.InsertText(f->getBeginLoc(), ", ", true, true);
                }
            }
            mRewriter.InsertText(f->getBeginLoc(), ">\n", true, true);
        }

        return true;
    }

    bool TraverseStmt(Stmt* e) {
        if (e && isa<Expr>(e) &&
                (cast<Expr>(e)->getType().getAsString() == "struct Sequence" ||
                 cast<Expr>(e)->getType().getAsString() == "struct Assign"))
        {
            Expr* expr = cast<Expr>(e);
            assert(expr);

            // find begin and end
            SourceLocation begin = expr->getBeginLoc();
            SourceLocation end = expr->getEndLoc();

            // find length of last token
            unsigned endOffset = Lexer::MeasureTokenLength(
                        end, mRewriter.getSourceMgr(), mRewriter.getLangOpts());
            endOffset++; // skip ';' after last token
            end = end.getLocWithOffset(endOffset);

            // process expression

            SequenceVisitor sv(mRewriter);
            sv.TraverseStmt(expr);

            if (sv.allVars().empty()) {
                // sometimes it can be empty (in headers)
                // but if it is empty in our code we should stop on error
                mRewriter.InsertText(begin, "!!! ``` !!!", true, true);
                return true;
            }

            mRewriter.InsertText(begin, "{\n", true, true);

            // put all rhs Vars into temporaries
            for (auto& e : sv.allVars()) {
                std::string ss = (Twine("Var ") + temporaryName(e) + " = " + e + ";\n").str();
                mRewriter.InsertText(begin, ss, true, true);
            }

            // find iteration space size
            auto sizeStr = (Twine("Size _size = ") + temporaryName(*sv.allVars().begin()) + ".range().size();\n").str();
            mRewriter.InsertText(begin, sizeStr, true, true);

            // add assertions
            for (auto& e : sv.allVars()) {
                std::string ss = (Twine("assert(") + temporaryName(e) + ".range().size() == _size);\n").str();
                mRewriter.InsertText(begin, ss, true, true);
            }

            // add expression lambda

            mRewriter.InsertText(begin, "auto _evaluate = [&](Int _i, Int _j, auto type) __attribute__((always_inline)) {\n", true, true);
            mRewriter.InsertText(begin, "using T = decltype(type);\n", true, true);

            // create temporaries for each loop iteration
            mRewriter.InsertText(begin, "// temporaries\n", true, true);
            for (auto& e : sv.temps) {
                std::string ss = (Twine("T ") + escapedName(e) + ";\n").str();
                mRewriter.InsertText(begin, ss, true, true);
            }
            
            mRewriter.InsertText(begin, "// variables\n", true, true);
            for (auto& e : sv.allVars()) {
                std::string ss = (Twine("T ") + escapedName(e) + ";\n").str();
                mRewriter.InsertText(begin, ss, true, true);
            }

            // load values from both lhs and rhs
            // we load lhs, because it is difficult to detect is it read only or read write
            mRewriter.InsertText(begin, "// loads\n", true, true);
            for (auto& e : sv.allVars()) {
                std::string ss = (Twine("loadFromPtr(") + escapedName(e) + ", " +
                                  "&" + temporaryName(e) + ".val(_i, _j));\n").str();
                mRewriter.InsertText(begin, ss, true, true);
            }

            // here the main expressions will be inserted
            // add end of line after it
            mRewriter.InsertText(end, "\n", true, true);
            
            // store values to lhs
            for (auto& e : sv.lhsVars) {
                std::string ss = (Twine("storeToPtr(") + 
                    "&" + temporaryName(e) + ".val(_i, _j)," +
                    escapedName(e) + ");\n").str();
                mRewriter.InsertText(end, ss, true, true);
            }
            
            // add closing brace of lambda
            mRewriter.InsertText(end, "};\n", true, true);

            // insert loops with call to lambda
            if (LoopsVectorization.getValue()) {
                mRewriter.InsertText(end, vecLoops, true, true);
            } else {
                mRewriter.InsertText(end, loops, true, true);
            }

            mRewriter.InsertText(end, "}\n", true, true);

            return true;
        } else {
            return RecursiveASTVisitor<MyASTVisitor>::TraverseStmt(e);
        }
    }

private:
    Rewriter& mRewriter;
};

class MyASTConsumer : public ASTConsumer {
public:
    MyASTConsumer(Rewriter& rewriter)
        : mVisitor(rewriter)
    {

    }

    bool HandleTopLevelDecl(DeclGroupRef dr) override {
        for (DeclGroupRef::iterator b = dr.begin(); b != dr.end(); ++b) {
            mVisitor.TraverseDecl(*b);
        }
        return true;
    }

private:
    MyASTVisitor mVisitor;
};

class MyFrontendAction : public ASTFrontendAction {
public:
    void EndSourceFileAction() override {
        SourceManager &sm = mRewriter.getSourceMgr();

        StringRef fname = sm.getFileEntryForID(sm.getMainFileID())->getName();
        llvm::errs() << "** EndSourceFileAction for: " << fname << "\n";

        std::error_code error_code;
        llvm::raw_fd_ostream outFile(MyToolOutput.getValue(), error_code,
                                     llvm::sys::fs::OF_None);
        mRewriter.getEditBuffer(sm.getMainFileID()).write(outFile);
    }

    std::unique_ptr<ASTConsumer> CreateASTConsumer(
            CompilerInstance& ci, StringRef file) override
    {
        llvm::errs() << "** Creating AST consumer for: " << file << "\n";
        mRewriter.setSourceMgr(ci.getSourceManager(), ci.getLangOpts());
        return llvm::make_unique<MyASTConsumer>(mRewriter);
    }
private:
    Rewriter mRewriter;
};

int main(int argc, const char **argv) {
    CommonOptionsParser op(argc, argv, MyToolCategory);
    ClangTool Tool(op.getCompilations(), op.getSourcePathList());

    return Tool.run(newFrontendActionFactory<MyFrontendAction>().get());
}

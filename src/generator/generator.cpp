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


const char loops[] = "for (Int _j = 0; _j < _size.y; _j++)\n"
                     "for (Int _i = 0; _i < _size.x; _i++)\n"
                     "_evaluate(_i, _j);\n";

//const char vecLoops[] =
//        "for (Int _j = 0; _j < _size.y; _j++) {\n"
//        "for (Int _i = 0; _i < _size.x / 4 * 4; _i+=4) {\n";

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

            std::string s = Lexer::getSourceText(
                        CharSourceRange::getTokenRange(e->getSourceRange()),
                        mRewriter.getSourceMgr(),
                        LangOptions());
            std::string es = escapedName(s);

            mRewriter.ReplaceText(e->getSourceRange(), es);

            mRewriter.InsertTextAfterToken(e->getEndLoc(), ".val(_i, _j)");

            // highly inefficient =)
            bool doNotAdd = false;
            for (auto& e : vars) {
                if (e.first == s || e.second == es) {
                    assert(e.first == s && e.second == es);
                    doNotAdd = true;
                }
            }

            if (!doNotAdd) {
                vars.push_back(std::make_pair(s, es));
            }

            return true;
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

    std::vector<std::pair<std::string, std::string>> vars;

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
            auto r = v->getTypeSourceInfo()->getTypeLoc().getSourceRange();
            mRewriter.ReplaceText(r, "Real");
        }
        return true;
    }

    bool VisitFunctionDecl(FunctionDecl *f) {
        if (f->getReturnType().getAsString() == "struct Expr") {
            mRewriter.ReplaceText(f->getReturnTypeSourceRange(), "Real");
        }
        for (unsigned i = 0; i < f->getNumParams(); i++) {
            ParmVarDecl* p = f->getParamDecl(i);
            if (p->getType().getAsString() == "struct Expr") {
                mRewriter.ReplaceText(p->getTypeSourceInfo()->getTypeLoc().getSourceRange(), "Real");
            }
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

            if (sv.vars.empty()) {
                // sometimes it can be empty (in headers)
                // but if it is empty in our code we should stop on error
                mRewriter.InsertText(begin, "!!! ``` !!!", true, true);
                return true;
            }

            mRewriter.InsertText(begin, "{\n", true, true);

            // put all Vars into temporaries
            for (auto& e : sv.vars) {
                std::string ss = (Twine("Var ") + (e.second) + " = " + e.first + ";\n").str();
                mRewriter.InsertText(begin, ss, true, true);
            }

            auto sizeStr = (Twine("Size _size = ") + sv.vars[0].second + ".range().size();\n").str();
            mRewriter.InsertText(begin, sizeStr, true, true);

            // add assertions
            for (auto& e : sv.vars) {
                std::string ss = (Twine("assert(") + e.second + ".range().size() == _size);\n").str();
                mRewriter.InsertText(begin, ss, true, true);
            }

            // add expression lambda

            mRewriter.InsertText(begin, "auto _evaluate = [&](Int _i, Int _j) [[gnu::always_inline]] {\n", true, true);
            mRewriter.InsertText(end, "\n};\n", true, true);

            // insert loops
            mRewriter.InsertText(end, loops, true, true);

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

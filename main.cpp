#include <stdio.h>
#include <string.h>
#include <boost/shared_ptr.hpp>
#include <Ptexture.h>
#include <PtexUtils.h>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>

#define PTEX_META_VERT_POSITIONS     "PtexVertPositions"
#define PTEX_META_FACE_VERT_COUNTS   "PtexFaceVertCounts"
#define PTEX_META_FACE_VERT_INDICES  "PtexFaceVertIndices"

class PtexFile {
public:
    PtexFile();
    ~PtexFile();
    bool Load(const char *filename);
    void ConcatMeta(std::vector<float> &points,
                    std::vector<int> &vertCounts,
                    std::vector<int> &vertIndices,
                    int vertOffset);
    void Write(PtexWriter *writer, int faceOffset,
               Ptex::DataType dataType, int numChannels, int reduce);
    int GetNumFaces() const { return _numFaces; }
    Ptex::DataType GetDataType() const { return _ptex->dataType(); }
    int GetNumChannels() const { return _ptex->numChannels(); }

private:
    PtexFile(const PtexFile &) {}
    PtexFile & operator = (const PtexFile &);

    PtexTexture *_ptex;
    int _numFaces;
    std::vector<float> _points;
    std::vector<int> _vertCounts;
    std::vector<int> _vertIndices;
};

typedef boost::shared_ptr<PtexFile> PtexFilePtr;
typedef std::vector<PtexFilePtr> PtexFileList;

PtexFile::PtexFile() : _ptex(NULL) {
}

PtexFile::~PtexFile() {
    if (_ptex) _ptex->release();
}

bool
PtexFile::Load(const char *filename) {
    Ptex::String ptexError;
    _ptex = PtexTexture::open(filename, ptexError, true);

    if (_ptex == NULL) {
        printf("Error: %s\n", ptexError.c_str());
        return false;
    }

    PtexMetaData* meta = _ptex->getMetaData();
    if (meta->numKeys() < 3) return false;

    const float* vp;
    const int *vi, *vc;
    int nvp, nvi, nvc;

    meta->getValue(PTEX_META_FACE_VERT_COUNTS, vc, nvc);
    if (nvc == 0)
        return false;

    meta->getValue(PTEX_META_VERT_POSITIONS, vp, nvp);
    if (nvp == 0)
        return false;

    meta->getValue(PTEX_META_FACE_VERT_INDICES, vi, nvi);
    if (nvi == 0)
        return false;

    for (int i = 0; i < nvp/3; ++i) {
        for (int j = 0; j < 3; ++j) {
            float v = vp[i*3+j];
            _points.push_back(v);
        }
    }
    const int *fv = vi;
    for (int i = 0, ptxidx = 0; i < nvc; ++i) {
        int nv = vc[i];

        _vertCounts.push_back(nv);
        for (int j = 0; j < nv; ++j) {
            _vertIndices.push_back(fv[j]);
        }
        fv += nv;
    }
    _numFaces = _ptex->numFaces();

    meta->release();
    return true;
}

void
PtexFile::ConcatMeta(std::vector<float> &points,
                     std::vector<int> &vertCounts,
                     std::vector<int> &vertIndices,
                     int vertOffset) {
    points.insert(points.end(), _points.begin(), _points.end());
    vertCounts.insert(vertCounts.end(),
                      _vertCounts.begin(), _vertCounts.end());

    for (std::vector<int>::iterator it = _vertIndices.begin();
         it != _vertIndices.end(); ++it) {
        vertIndices.push_back(*it + vertOffset);
    }
}

void
PtexFile::Write(PtexWriter *writer, int faceOffset, Ptex::DataType dataType,
                int numChannels, int reduce) {
    int srcBpp = Ptex::DataSize(_ptex->dataType()) * _ptex->numChannels();
    int dstBpp = Ptex::DataSize(dataType) * numChannels;

    float *buffer = new float[numChannels];


    for (int i = 0; i < _numFaces; ++i) {
        Ptex::FaceInfo info = _ptex->getFaceInfo(i);
        // fix adjacency face index
        int adjfaces[4];
        for (int k = 0; k < 4; ++k) {
            adjfaces[k] = info.adjface(k);
            if (adjfaces[k] != -1) adjfaces[k] += faceOffset;
        }
        info.setadjfaces(adjfaces[0], adjfaces[1], adjfaces[2], adjfaces[3]);

        if (reduce != 0) {
            info.res = Ptex::Res(std::max(info.res.ulog2 - reduce, 0),
                                 std::max(info.res.vlog2 - reduce, 0));
            int dstSize = dstBpp * info.res.size();
            unsigned char *dstBuffer = new unsigned char[dstSize];
            int index = 0;
            for (int v = 0; v < info.res.v(); ++v) {
                for (int u = 0; u < info.res.u(); ++u) {
                    _ptex->getPixel(i, u, v, buffer, 0, numChannels, info.res);

                    Ptex::ConvertFromFloat(dstBuffer + dstBpp * index, buffer,
                                           dataType, numChannels);
                    index++;
                }
            }
            writer->writeFace(faceOffset + i, info, dstBuffer);
            delete[] dstBuffer;
        } else {
            int numTexels = info.res.size();

            // get face texels
            int srcSize = srcBpp * numTexels;
            unsigned char *srcBuffer = new unsigned char[srcSize];
            _ptex->getData(i, srcBuffer, 0);

            if (dataType == _ptex->dataType() &&
                numChannels == _ptex->numChannels()) {
                // write face directly
                writer->writeFace(faceOffset + i, info, srcBuffer);
            } else {
                // needs format conversion
                int dstSize = dstBpp * numTexels;
                unsigned char *dstBuffer = new unsigned char[dstSize];

                float color[4];
                for (int j = 0; j < numTexels; ++j) {
                    Ptex::ConvertToFloat(color, srcBuffer + srcBpp * j,
                                         _ptex->dataType(),
                                         _ptex->numChannels());
                    Ptex::ConvertFromFloat(dstBuffer + dstBpp * j, color,
                                           dataType, numChannels);
                }
                // write face
                writer->writeFace(faceOffset + i, info, dstBuffer);

                delete[] dstBuffer;
            }
            delete[] srcBuffer;
        }
    }
    delete[] buffer;
}

// -----------------------------------------------------------------------------

static void usage(const char *program) {
    printf("Usage: %s [options] -o <outfile.ptx> "
           "<infile1.ptx> <infile2.ptx> ...\n", program);
    printf("Options: -d uint8|uint16|half|float    out data type\n");
    printf("         -n 1, 2, 3, 4                 out num channels\n");
    printf("         -e <out.obj>                  export merged obj file\n");
    printf("         -r reduce_level               reduce resolution "
           "(1=half, 2=quater, ..)\n");
}

int main(int argc, char *argv[]) {
    std::string outfile;
    std::string objfile;
    std::vector<std::string> infiles;
    const char *dataTypeStr[] = { "uint8", "uint16", "half", "float" };
    int dataType = -1;
    int numChannels = -1;
    int reduce = 0;

    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
            case 'o':
                outfile = argv[++i];
                break;
            case 'd':
            {
                const char *dt = argv[++i];
                for (int j = 0; j < 4; ++j) {
                    if (!strcmp(dataTypeStr[j], dt)) dataType = j;
                }
                if (dataType == -1) {
                    usage(argv[0]);
                    return 1;
                }
                break;
            }
            case 'n':
                numChannels = atoi(argv[++i]);
                break;
            case 'e':
                objfile = argv[++i];
                break;
            case 'r':
                reduce = atoi(argv[++i]);
                break;
            default:
                usage(argv[0]);
                return 1;
            }
        } else {
            infiles.push_back(argv[i]);
        }
    }

    if (outfile.empty()) {
        usage(argv[0]);
        return 1;
    }

    PtexFileList ptexlist;
    std::vector<float> points;
    std::vector<int> vertCounts, vertIndices;
    int faceBase = 0;
    int numTotalFaces = 0;

    // read ptx files
    for (int i = 0; i < infiles.size(); ++i) {
        PtexFilePtr ptex(new PtexFile());
        if (ptex->Load(infiles[i].c_str()) == false) {
            printf("Error: reading %s\n", infiles[i].c_str());
            return 1;
        }
        int vertOffset = points.size()/3;
        ptex->ConcatMeta(points, vertCounts, vertIndices, vertOffset);
        numTotalFaces += ptex->GetNumFaces();

        if (dataType == -1)
            dataType = ptex->GetDataType();
        if (numChannels == -1)
            numChannels = ptex->GetNumChannels();

        if (dataType != ptex->GetDataType() ||
            numChannels != ptex->GetNumChannels()) {
            printf("Warning: converting %s from %s/%d to %s/%d\n",
                   infiles[i].c_str(),
                   dataTypeStr[ptex->GetDataType()],
                   ptex->GetNumChannels(),
                   dataTypeStr[dataType],
                   numChannels);
        }

        ptexlist.push_back(ptex);
    }

    // merged obj export
    if (not objfile.empty()) {
        FILE *fp = fopen(objfile.c_str(), "w");
        size_t np = points.size()/3;
        for (size_t i = 0; i < np; ++i) {
            fprintf(fp, "v %f %f %f\n",
                    points[i*3+0], points[i*3+1], points[i*3+2]);
        }
        size_t nf = vertCounts.size();
        int index = 0;
        for (size_t i = 0; i < nf; ++i) {
            int nv = vertCounts[i];
            fprintf(fp, "f");
            for (int j = 0; j < nv; ++j) {
                fprintf(fp, " %d", vertIndices[index++]+1);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    // merged ptex write
    Ptex::String ptexError;
    int alphachan = numChannels == 4 ? 3 : -1;
    PtexWriter *writer = PtexWriter::open(outfile.c_str(), Ptex::mt_quad,
                                          (Ptex::DataType)dataType, numChannels,
                                          alphachan, numTotalFaces, ptexError,
                                          true);

    writer->writeMeta(PTEX_META_VERT_POSITIONS,
                      &points[0], points.size());
    writer->writeMeta(PTEX_META_FACE_VERT_COUNTS,
                      &vertCounts[0], vertCounts.size());
    writer->writeMeta(PTEX_META_FACE_VERT_INDICES,
                      &vertIndices[0], vertIndices.size());

    int faceOffset = 0;
    for (PtexFileList::iterator it = ptexlist.begin();
         it != ptexlist.end(); ++it) {
        (*it)->Write(writer, faceOffset,
                     (Ptex::DataType)dataType, numChannels, reduce);
        faceOffset += (*it)->GetNumFaces();
    }
    if (writer->close(ptexError) == false) {
        printf("Error: writing %s\n", outfile.c_str());
        printf("%s\n", ptexError.c_str());
        return 1;
    }

    writer->release();

    return 0;
}

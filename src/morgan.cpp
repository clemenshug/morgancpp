#include <Rcpp.h>
#include <zlib.h>
#include "lz4/lz4hc.h"
#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

using Fingerprint = std::array<std::uint64_t, 32>;
using FingerprintName = std::int32_t;
using FingerprintN = std::uint64_t;

enum {
  MESSAGE_MAX_BYTES   = 1024 * 64
};

#define CMPBUFSIZE (LZ4_COMPRESSBOUND(MESSAGE_MAX_BYTES))
#define DECODER_RING_BUFFER_SIZE (LZ4_DECODER_RING_BUFFER_SIZE(MESSAGE_MAX_BYTES))

const int ZERO = 0;

int lz4_compress(char* out_buffer, const char* in_buffer, std::size_t in_size, const int compression_level) {
  LZ4_streamHC_t* lz4Stream = LZ4_createStreamHC();
  LZ4_resetStreamHC_fast(lz4Stream, compression_level);

  std::size_t inpOffset = 0;
  std::size_t outOffset = 0;
  std::size_t cmpBufOffset = 0;

  for(;;) {
    char cmpBuf[CMPBUFSIZE];;

    int inpBytes = in_size - inpOffset;
    if (0 == inpBytes)
      break;
    if (inpBytes > MESSAGE_MAX_BYTES)
      inpBytes = MESSAGE_MAX_BYTES;

    const int cmpBytes = LZ4_compress_HC_continue(
      lz4Stream,
      in_buffer + inpOffset,
      cmpBuf,
      inpBytes,
      CMPBUFSIZE
    );
    if(cmpBytes <= 0)
      break;
    inpOffset += inpBytes;

    memcpy(out_buffer + outOffset, reinterpret_cast<const char *>(&cmpBytes), sizeof(int));
    outOffset += sizeof(int);

    memcpy(out_buffer + outOffset, cmpBuf, cmpBytes);
    outOffset += cmpBytes;
  }

  memcpy(out_buffer + outOffset, reinterpret_cast<const char *>(&ZERO), sizeof(int));
  outOffset += sizeof(int);

  return outOffset;
}


int lz4_decompress(char* out_buffer, const char* in_buffer)
{
  LZ4_streamDecode_t* lz4Stream = LZ4_createStreamDecode();

  std::size_t inpOffset = 0;
  std::size_t outOffset = 0;
  static char decBuf[DECODER_RING_BUFFER_SIZE];
  std::size_t decBufOffset = 0;

  // Rcpp::Rcout << "Using a " << DECODER_RING_BUFFER_SIZE << " byte decoding ring buffer\n";

  for(;;) {
    int const inpBytes = *reinterpret_cast<const int *>(in_buffer + inpOffset);
    inpOffset += sizeof(int);
    // Rcpp::Rcout << "Reading " << inpBytes << " bytes\n";
    if (inpBytes <= 0)
      break;

    char* const decPtr = &decBuf[decBufOffset];
    // Rcpp::Rcout << "Decode buffer offset is " << decBufOffset << " bytes\n";
    // Rcpp::Rcout << "Input buffer offset is " << inpOffset << " bytes\n";
    // Rcpp::Rcout << "Output buffer offset is " << outOffset << " bytes\n";
    int const decBytes = LZ4_decompress_safe_continue(
      lz4Stream,
      in_buffer + inpOffset,
      decPtr,
      inpBytes,
      MESSAGE_MAX_BYTES
    );
    // Rcpp::Rcout << "Decompressed into " << decBytes << " bytes\n";
    if(decBytes <= 0) {
      Rcpp::Rcout << "Decompression error. inpOffset: " << inpOffset << "\noutOffset: " << outOffset
                  << "\ndecBufOffset: " << decBufOffset << "\ninpBytes: " << inpBytes
                  << "\ndecBytes: " << decBytes << "\nDECODER_RING_BUFFER_SIZE: " << DECODER_RING_BUFFER_SIZE
                  << "\nMESSAGE_MAX_BYTES: " << MESSAGE_MAX_BYTES << "\n";
      break;
    }


    inpOffset += inpBytes;

    memcpy(out_buffer + outOffset, decPtr, decBytes);
    outOffset += decBytes;

    decBufOffset += decBytes;
    // Wraparound the ringbuffer offset
    if(decBufOffset >= DECODER_RING_BUFFER_SIZE - MESSAGE_MAX_BYTES) {
      // Rcpp::Rcout << "Wrap ring buffer\n";
      decBufOffset = 0;
    }
  }

  return outOffset;
}

// Parse a single hexadecimal character to its integer representation.
int parse_hex_char(const char& c) {
  int v;
  if ((c >= '0') && (c <= '9')) {
    v = (c - '0');
  } else if ((c >= 'A') && (c <= 'F')) {
    v = (c - 'A' + 10);
  } else if ((c >= 'a') && (c <= 'f')) {
    v = (c - 'a' + 10);
  } else {
    ::Rf_error("Hex string may only contain characters in [0-9A-F]");
  }
  return v;
}

// Convert raw byte string to fingerprint.
Fingerprint raw2fp(const std::string& raw) {
  if (raw.length() != 256) {
    ::Rf_error("Input raw string must be of length 256");
  }
  Fingerprint fp;
  std::memcpy(fp.data(), &raw[0], sizeof(fp));
  return fp;
}

// Convert ASCII hex string to fingerprint.
Fingerprint hex2fp(const std::string& hex) {
  if (hex.length() != 512) {
    ::Rf_error("Input hex string must be of length 512");
  }
  // Convert hex to raw bytes
  std::string raw(256, '\0');
  auto hi = std::cbegin(hex);
  for (auto& ri : raw) {
    ri = parse_hex_char(*hi++) | (parse_hex_char(*hi++) << 4);
  }
  return raw2fp(raw);
}


// Compute Jaccard similarity of two fingerprints
double jaccard_fp(const Fingerprint& f1, const Fingerprint& f2) {
  int count_and = 0, count_or = 0;
  auto i1 = std::cbegin(f1), i2 = std::cbegin(f2);
  while (i1 != std::cend(f1)) {
    count_and += __builtin_popcountll(*i1 & *i2);
    count_or += __builtin_popcountll(*i1 | *i2);
    ++i1;
    ++i2;
  }
  return static_cast<double>(count_and) / count_or;
}


//' Tanimoto similarity between two Morgan fingerprints
//'
//' Computes Tanimoto similarity between two hexadecimal strings
//'
//' @param s1 Hexadecimal string of length 512
//' @param s2 Hexadecimal string of length 512
//' @return Jaccard similarity over the bits representing individual keys
//' @export
// [[Rcpp::export]]
double tanimoto(const std::string& s1, const std::string& s2) {
  const Fingerprint& fp1 = hex2fp(s1);
  const Fingerprint& fp2 = hex2fp(s2);
  return jaccard_fp(fp1, fp2);
}


//' @name MorganFPS
//' @title Morgan fingerprints collection
//' @description Efficient structure for storing a set of Morgan fingerprints
//' @field new Constructor. Accepts either a vector of fingerprints in hexadecimal
//'   format or a path to a binary file of fingerprints using the argument
//'   `from_file = TRUE`
//' @field tanimoto (i,j) similarity between fingerprints i and j
//' @field tanimoto_all (i) similarity between fingerprint i and all others
//' @field tanimoto_ext (s) similarity between external hexadecimal string s and all
//'    fingerprints in the collection
//' @field save_file (path, compression_level) Save fingerprints to file in binary format
//' @field size number of bytes used to store the fingerprints
//' @importFrom Rcpp cpp_object_initializer
//' @export
class MorganFPS {

public:

  // Constructor accepts a named character vector of hex strings
  MorganFPS(Rcpp::CharacterVector fps_hex) {
    size_t n = fps_hex.length();
    Rcpp::RObject passed_names = fps_hex.names();
    fps.reserve(n);
    fp_names.reserve(n);
    if(passed_names.isNULL()) {
      for (FingerprintName i = 0; i < n; i++) {
        fp_names.push_back(i + 1);
      }
    } else {
      Rcpp::CharacterVector passed_names_vec = Rcpp::as<Rcpp::CharacterVector>(passed_names);
      if (Rcpp::unique(passed_names_vec).length() != n)
        Rcpp::stop("Names must be unique");
      for (Rcpp::CharacterVector::iterator i = passed_names_vec.begin(); i != passed_names_vec.end(); i++) {
        fp_names.push_back(std::stoll(std::string(*i)));
      }
      Rcpp::IntegerVector idx = Rcpp::seq_along(passed_names_vec) - 1;
      std::sort(idx.begin(), idx.end(), [&](int i, int j){return fp_names[i] < fp_names[j];});
      std::vector<FingerprintName> fp_names_sorted;
      fp_names_sorted.reserve(n);
      for (auto &i : idx) {
        fp_names_sorted.push_back(fp_names[i]);
      }
      fp_names = fp_names_sorted;
      fps_hex = fps_hex[idx];
    }
    for (Rcpp::CharacterVector::iterator i = fps_hex.begin(); i != fps_hex.end(); i++) {
      fps.push_back(hex2fp(std::string(*i)));
    }
  }

  // Constructor accepts a file path to load fingerprints from binary file
  MorganFPS(const std::string& filename, const bool from_file) {
    std::ifstream in_stream;
    in_stream.open(filename, std::ios::in | std::ios::binary | std::ios::ate);
    std::streampos file_size = in_stream.tellg();
    std::vector<char> compressed_buffer;

    compressed_buffer.resize(9);
    in_stream.seekg (0, in_stream.beg);
    in_stream.read(compressed_buffer.data(), 9);
    std::string magic (compressed_buffer.begin(), compressed_buffer.end());
    if (magic != "MORGANFPS") {
      Rcpp::stop("File is incompatible, doesn't start with 'MORGANFPS': '%s'", magic.c_str());
    }

    FingerprintN n;
    in_stream.read(reinterpret_cast<char*>(&n), sizeof(FingerprintN));
    Rcpp::Rcout << "Reading " << n << " fingerprints from file\n";
    fps.resize(n);
    fp_names.resize(n);

    int size_next_block;
    in_stream.read(reinterpret_cast<char*>(&size_next_block), sizeof(int));
    Rcpp::Rcout << "Fingerprint block has " << size_next_block << " bytes\n";

    compressed_buffer.resize(size_next_block);
    in_stream.read(compressed_buffer.data(), size_next_block);
    Rcpp::Rcout << "Compressed fingerprints read\n";
    int bytes_decompressed = lz4_decompress(
      reinterpret_cast<char*>(fps.data()),
      compressed_buffer.data()
    );
    if (bytes_decompressed != n * sizeof(Fingerprint)) {
      Rcpp::stop(
        "Decompression error in fingerprints:\nExpected bytes: %i\nReceived bytes: %i",
        n * sizeof(Fingerprint),
        bytes_decompressed
      );
    }
    Rcpp::Rcout << "Fingerprints decompressed\n";

    in_stream.read(reinterpret_cast<char*>(&size_next_block), sizeof(int));
    Rcpp::Rcout << "Names block has " << size_next_block << " bytes\n";

    compressed_buffer.resize(size_next_block);
    in_stream.read(compressed_buffer.data(), size_next_block);
    bytes_decompressed = lz4_decompress(
      reinterpret_cast<char*>(fp_names.data()),
      compressed_buffer.data()
    );
    if (bytes_decompressed != n * sizeof(FingerprintName)) {
      Rcpp::stop(
        "Decompression error in names:\nExpected bytes: %i\nReceived bytes: %i",
        n * sizeof(Fingerprint),
        bytes_decompressed
      );
    }
    Rcpp::Rcout << "Names decompressed\n";
  }

  // Tanimoto similarity between drugs i and j
  double tanimoto(FingerprintName &i, FingerprintName &j) {
    return jaccard_fp(fp_index(i), fp_index(j));
  }

  // Tanimoto similarity of drug i to every other drug
  Rcpp::DataFrame tanimoto_all(FingerprintName &i) {
    const Fingerprint& fp_other = fp_index(i);
    Rcpp::NumericVector res(fps.size());
    for (int i = 0; i < fps.size(); i++) {
      res[i] = jaccard_fp(fps[i], fp_other);
    }
    return Rcpp::DataFrame::create(
      Rcpp::Named("id") = fp_names,
      Rcpp::Named("structural_similarity") = res
    );
  }

  // Tanimoto similarity of an external drug to every other drug
  //   in the collection
  Rcpp::DataFrame tanimoto_ext(const std::string& other) {
    const Fingerprint& fp_other = hex2fp(other);
    Rcpp::NumericVector res(fps.size());
    for (int i = 0; i < fps.size(); i++) {
      res[i] = jaccard_fp(fps[i], fp_other);
    }
    return Rcpp::DataFrame::create(
      Rcpp::Named("id") = fp_names,
      Rcpp::Named("structural_similarity") = res
    );
  }

  void save_file(const std::string& filename) {
    save_file(filename, 8);
  }

  // Save binary fp file
  void save_file(const std::string& filename, const int& compression_level=50) {
    if (compression_level < 0 || compression_level > 100)
      Rcpp::stop("Compression level must be between 0 and 100");

    FingerprintN n = fps.size();
    Rcpp::Rcout << "Wrinting " << n << " fingerprints\n";

    std::ofstream out_stream;
    out_stream.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);

    std::vector<char> out_buffer;

    out_stream.write("MORGANFPS", 9);
    out_stream.write(reinterpret_cast<char*>(&n), sizeof(FingerprintN));

    out_buffer.resize(LZ4_COMPRESSBOUND(fps.size() * sizeof(Fingerprint)));
    const int fingerprints_compressed = lz4_compress(
      out_buffer.data(),
      reinterpret_cast<char *>(fps.data()),
      fps.size() * sizeof(Fingerprint),
      compression_level
    );

    Rcpp::Rcout << "Fingerprints compressed " << fingerprints_compressed << " bytes\n";
    // Save number of bytes of the compressed data. Important for finding
    // second block with names for decompression
    out_stream.write(reinterpret_cast<const char*>(&fingerprints_compressed), sizeof(int));
    out_stream.write(out_buffer.data(), fingerprints_compressed);
    Rcpp::Rcout << "Wrote fingerprints\n";

    out_buffer.resize(LZ4_COMPRESSBOUND(fp_names.size() * sizeof(FingerprintName)));
    const int names_compressed = lz4_compress(
      out_buffer.data(),
      reinterpret_cast<char *>(fp_names.data()),
      fp_names.size() * sizeof(FingerprintName),
      compression_level
    );

    Rcpp::Rcout << "Names compressed " << names_compressed << " bytes\n";
    // Save number of bytes of the compressed data. Important for finding
    // second block with names for decompression
    out_stream.write(reinterpret_cast<const char*>(&names_compressed), sizeof(int));
    out_stream.write(out_buffer.data(), names_compressed);
    Rcpp::Rcout << "Wrote Names\n";

    out_stream.close();
  }

  // Size of the dataset
  int size() {
    return fps.size() * sizeof(Fingerprint);
  }

  std::vector<Fingerprint> fps;
  std::vector<FingerprintName> fp_names;

private:

  Fingerprint& fp_index(FingerprintName &x) {
    auto fp_pt = std::lower_bound(fp_names.begin(), fp_names.end(), x);
    if (*fp_pt != x)
      Rcpp::stop("Fingerprint %s not found", std::to_string(x).c_str());
    return fps.at(fp_pt - fp_names.begin());
  }

};

// Expose all relevant classes through an Rcpp module
RCPP_EXPOSED_CLASS(MorganFPS)
RCPP_MODULE(morgan_cpp) {

  using namespace Rcpp;

  class_<MorganFPS>( "MorganFPS" )
    .constructor< CharacterVector >("Construct fingerprint collection from vector of fingerprints")
    .constructor< std::string, bool >("Construct fingerprint collection from binary file")
    .method("size", &MorganFPS::size, "Size of the data in bytes")
    .method("tanimoto", &MorganFPS::tanimoto,
	    "Similarity between two fingerprints in the collection")
    .method("tanimoto_all", &MorganFPS::tanimoto_all,
	    "Similarity of a fingerprint against all other fingerprints in the collection")
    .method("tanimoto_ext", &MorganFPS::tanimoto_ext,
	    "Similarity of an external fingerprints against all fingerprints in the collection")
    .method("save_file", (void (MorganFPS::*)(const std::string&, const int&)) (&MorganFPS::save_file),
	    "Save fingerprints to file in binary format")
    .method("save_file", (void (MorganFPS::*)(const std::string&)) (&MorganFPS::save_file),
	    "Save fingerprints to file in binary format")
    .field_readonly("fingerprints", &MorganFPS::fps)
    .field_readonly("names", &MorganFPS::fp_names)
    ;
}

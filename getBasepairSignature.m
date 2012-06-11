function [signature] = getBasepairSignature(filename)

    load(filename);
    signature = Search.Signature;

end
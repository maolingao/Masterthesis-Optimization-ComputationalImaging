function dataStr_extended = extendDataStr(dataStr,len)

    dataStr_extended.errs = extend2len(dataStr.errs,len);
    dataStr_extended.rerrs = extend2len(dataStr.rerrs,len);
    dataStr_extended.time = extend2len(dataStr.time,len);
end
def pipe_dec(original):
    class PipeInto(object):
        data = {'function': original}

        def __init__(self, *args, **kwargs):
            self.data['args'] = args
            self.data['kwargs'] = kwargs

        def __rrshift__(self, other):
            return self.data['function'](
                self,
                other,
                *self.data['args'],
                **self.data['kwargs']
            )

    return PipeInto
